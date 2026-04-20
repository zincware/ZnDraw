import { ErrorOutline as ErrorIcon } from "@mui/icons-material";
import { Box, CircularProgress, Typography } from "@mui/material";
import { useColorScheme } from "@mui/material/styles";
import type { IDockviewPanelProps } from "dockview-react";
// @ts-expect-error - plotly.js internal module lacks type declarations
import { decodeTypedArraySpec } from "plotly.js/src/lib/array.js";
import type * as PlotlyJS from "plotly.js-dist-min";
import Plotly from "plotly.js-dist-min";
import { useCallback, useEffect, useMemo, useRef } from "react";
import { useFigure } from "../hooks/useFigures";
import { useStepControl } from "../hooks/useStepControl";
import { useAppStore } from "../store";

// Extend Plotly types to include custom meta property for step interactions
declare module "plotly.js-dist-min" {
	interface PlotData {
		meta?: {
			interactions?: Array<{
				click?: string;
				select?: string;
				hover?: string;
			}>;
		};
	}
}

/**
 * Tracks marker position for each interaction (step, geometry selection, etc.)
 * Used to efficiently update marker positions without full re-renders.
 */
interface MarkerTrackInfo {
	dataTraceIndex: number; // Index of original trace in plotData
	markerTraceIndex: number; // Index of marker trace in plotData
	xData: number[]; // x values of the original trace
	yData: number[]; // y values of the original trace
	customdata: any[]; // customdata array for finding points at frame/selection
	interactionType: "step" | "selection"; // Type of interaction
	interactionKey: string; // "step" or geometry name (e.g., "particles", "forces")
	customdataIndex: number; // Which dimension in customdata array to use
}

/**
 * Marker styling configuration for different interaction types.
 */
const MARKER_CONFIGS = {
	step: {
		size: 12,
		color: "red",
		borderColor: "darkred",
		borderWidth: 2,
	},
	selection: {
		size: 8,
		color: "rgb(100, 150, 200)", // Lighter steel blue for selections
		borderColor: "rgb(50, 100, 160)",
		borderWidth: 1.5,
	},
} as const;

/**
 * Converts typed array specs (bdata + dtype) to regular arrays.
 * Handles Plotly's binary data encoding and reshaping for multi-dimensional data (arbitrary shapes).
 */
function convertToArray(data: any): any[] {
	if (!data) {
		return [];
	}

	if (Array.isArray(data)) {
		return data;
	}

	if (typeof data === "object" && "bdata" in data && "dtype" in data) {
		try {
			const flatArray = decodeTypedArraySpec(data);

			if (!Array.isArray(flatArray)) {
				// Convert TypedArray to regular array
				if (
					flatArray &&
					typeof flatArray === "object" &&
					flatArray.length !== undefined
				) {
					return Array.from(flatArray);
				}
				return [];
			}

			// If data has a shape parameter, reshape the array to match dimensions
			if (data.shape) {
				const shapeStr = data.shape;
				// Parse shape like "67, 2" to get all dimensions
				const dims = shapeStr
					.split(",")
					.map((s: string) => Number.parseInt(s.trim()));

				if (dims.length >= 1) {
					// Recursively reshape for arbitrary dimensions
					const result = reshapeFlat(flatArray, dims, 0);
					return result;
				}
			}

			return flatArray;
		} catch (_e) {
			return [];
		}
	}

	return [];
}

/**
 * Recursively reshape a flat array into an n-dimensional array.
 *
 * Parameters
 * ----------
 * flatArray : any[]
 *     The flat array to reshape.
 * dims : number[]
 *     Array of dimensions [d0, d1, d2, ...].
 * dimIndex : number
 *     Current dimension being processed.
 */
function reshapeFlat(flatArray: any[], dims: number[], dimIndex: number): any {
	if (dimIndex === dims.length - 1) {
		// Last dimension - return slice of flatArray
		const size = dims[dimIndex];
		return flatArray.slice(0, size);
	}

	// Recursive case - create array of sub-arrays
	const size = dims[dimIndex];
	const subArraySize = dims.slice(dimIndex + 1).reduce((a, b) => a * b, 1);
	const result: any[] = [];

	for (let i = 0; i < size; i++) {
		const start = i * subArraySize;
		const end = start + subArraySize;
		const subArray = reshapeFlat(
			flatArray.slice(start, end),
			dims,
			dimIndex + 1,
		);
		result.push(subArray);
	}

	return result;
}

/**
 * Extracts a single numeric value from a customdata entry at a specific dimension.
 * This handles various complex data structures that Plotly might return.
 */
function getCustomDataValue(
	customDataPoint: any,
	dimensionIndex: number,
): number | undefined {
	// Skip if customdata point is null/undefined
	if (customDataPoint === null || customDataPoint === undefined) {
		return undefined;
	}

	let pointValue: any; // Use 'any' to handle intermediate values

	// Handle array-like customdata (most common case for (m,n) shapes)
	if (Array.isArray(customDataPoint)) {
		const element = customDataPoint[dimensionIndex];

		// Skip if dimension doesn't exist or is null/undefined
		if (element === null || element === undefined) {
			return undefined;
		}

		// Extract value from element
		if (typeof element === "number") {
			pointValue = element;
		} else if (Array.isArray(element)) {
			// Nested array - get first element
			pointValue = element[0];
		} else if (
			typeof element === "object" &&
			element.length !== undefined &&
			"buffer" in element
		) {
			// TypedArray - get first value
			pointValue = element[0];
		} else if (typeof element === "object" && element !== null) {
			// Object with numeric keys like {"0": 8, "1": 8}
			const numKeys = Object.keys(element).filter(
				(k) => !Number.isNaN(Number(k)),
			);
			if (numKeys.length > 0) {
				pointValue = (element as Record<string, any>)[numKeys[0]];
			}
		}
	} else if (typeof customDataPoint === "number") {
		// Scalar value - only match for dimension 0
		if (dimensionIndex === 0) {
			pointValue = customDataPoint;
		}
	} else if (typeof customDataPoint === "object") {
		// Object with numeric keys - handle for arbitrary dimension
		const numKeys = Object.keys(customDataPoint).filter(
			(k) => !Number.isNaN(Number(k)),
		);
		if (numKeys.length > 0) {
			const sortedKeys = numKeys.sort((a, b) => Number(a) - Number(b));
			const keyAtDim = sortedKeys[dimensionIndex];
			if (keyAtDim !== undefined) {
				pointValue = (customDataPoint as Record<string, any>)[keyAtDim];
			}
		}
	}

	// Ensure we return a number or undefined
	if (pointValue === null || pointValue === undefined) {
		return undefined;
	}

	const numericValue = Number(pointValue);
	return Number.isNaN(numericValue) ? undefined : numericValue;
}

/**
 * Find data points matching a specific value in a customdata dimension.
 * Returns array of {x, y} positions for all matching points.
 */
function findPointsWithValue(
	track: MarkerTrackInfo,
	value: number,
	dimensionIndex = 0,
): Array<{ x: number; y: number }> {
	const points: Array<{ x: number; y: number }> = [];

	// Case 1: customdata exists - use it to find matching values
	if (track.customdata.length > 0) {
		for (let i = 0; i < track.customdata.length; i++) {
			const customDataPoint = track.customdata[i];

			const pointValue = getCustomDataValue(customDataPoint, dimensionIndex);

			// biome-ignore lint/suspicious/noDoubleEquals: intentional loose equality — plotly values may be string or number, coercion is required
			if (pointValue != null && pointValue == value) {
				points.push({
					x: track.xData[i] ?? i,
					y: track.yData[i],
				});
			}
		}
	} else {
		// Case 2: no customdata - use index as value (1D plots)
		if (value >= 0 && value < track.yData.length) {
			points.push({
				x: track.xData[value] ?? value,
				y: track.yData[value],
			});
		}
	}

	return points;
}

/**
 * Find a data point at a specific frame in a trace (uses dimension 0 of customdata).
 * Returns the x, y position or null if no point exists at that frame.
 */
function findPointAtFrame(
	track: MarkerTrackInfo,
	frame: number,
): { x: number; y: number } | null {
	const points = findPointsWithValue(track, frame, track.customdataIndex);
	return points.length > 0 ? points[0] : null;
}

interface PlotViewParams {
	figureKey: string;
}

export function PlotView(props: IDockviewPanelProps<PlotViewParams>) {
	const figureKey = props.params?.figureKey ?? "";

	// ===== DOM REFERENCES =====
	// Plotly extends the HTML element at runtime with event methods
	const plotContainer = useRef<PlotlyJS.PlotlyHTMLElement | null>(null);

	// ===== MARKER TRACKING =====
	// Stores info about marker traces for efficient frame updates
	const markerTracks = useRef<MarkerTrackInfo[]>([]);

	// ===== EVENT LISTENER REFS =====
	// Keep references to handlers for proper cleanup
	const onClickRef = useRef<((data: PlotlyJS.PlotMouseEvent) => void) | null>(
		null,
	);
	const onSelectedRef = useRef<
		((data: PlotlyJS.PlotSelectionEvent) => void) | null
	>(null);
	const onDeselectRef = useRef<(() => void) | null>(null);

	// ===== DATA FETCHING =====
	const {
		data: figureResponse,
		isLoading,
		isError,
		error,
	} = useFigure(figureKey, { enabled: !!figureKey });

	// ===== STORE ACTIONS & STATE =====
	const updateSelectionForGeometry = useAppStore(
		(state) => state.updateSelectionForGeometry,
	);
	const updateFrameSelection = useAppStore(
		(state) => state.updateFrameSelection,
	);
	const currentFrame = useAppStore((state) => state.currentFrame);
	const frame_selection = useAppStore((state) => state.frame_selection);
	const selections = useAppStore((state) => state.selections);
	const { setStep } = useStepControl();
	const { mode } = useColorScheme();

	// ===== MEMOIZED: Parsed Plotly JSON =====
	const plotlyJson = useMemo(() => {
		if (
			figureResponse?.figure?.type === "plotly" &&
			figureResponse.figure.data
		) {
			try {
				return JSON.parse(figureResponse.figure.data);
			} catch (e) {
				console.error("Failed to parse Plotly JSON for figure:", {
					figureKey,
					error: e instanceof Error ? e.message : String(e),
				});
				return null;
			}
		}
		return null;
	}, [figureResponse?.figure, figureKey]);

	// ===== MEMOIZED: Layout with autosize and dark mode support =====
	const plotLayout = useMemo(() => {
		const isDark = mode === "dark";
		const darkModeOverrides = isDark
			? {
					paper_bgcolor: "rgba(30, 30, 30, 1)",
					plot_bgcolor: "rgba(30, 30, 30, 1)",
					font: { color: "#e0e0e0" },
					xaxis: {
						...plotlyJson?.layout?.xaxis,
						gridcolor: "rgba(100, 100, 100, 0.3)",
						zerolinecolor: "rgba(150, 150, 150, 0.5)",
					},
					yaxis: {
						...plotlyJson?.layout?.yaxis,
						gridcolor: "rgba(100, 100, 100, 0.3)",
						zerolinecolor: "rgba(150, 150, 150, 0.5)",
					},
				}
			: {};

		return {
			...plotlyJson?.layout,
			...darkModeOverrides,
			autosize: true,
		};
	}, [plotlyJson?.layout, mode]);

	// ===== STABLE CALLBACKS: Event Handlers =====

	const onPlotClick = useCallback(
		(data: PlotlyJS.PlotMouseEvent) => {
			const aggregatedSelections = new Map<string, number[]>();
			let stepFrame: number | undefined;

			data.points.forEach((point) => {
				if (
					point.data?.meta?.interactions &&
					Array.isArray(point.data.meta.interactions)
				) {
					// For each interaction dimension, extract the corresponding customdata value
					point.data.meta.interactions.forEach(
						(
							interaction:
								| { click?: string; select?: string; hover?: string }
								| null
								| undefined,
							dimensionIndex: number,
						) => {
							// Skip if interaction is null/undefined
							if (!interaction) return;

							// Skip if interaction has no click action
							if (!interaction.click) return;

							const numericValue = getCustomDataValue(
								point.customdata,
								dimensionIndex,
							);

							// Skip if customdata value is null/undefined
							if (numericValue === undefined) return;

							if (interaction.click === "step") {
								stepFrame = numericValue;
							} else if (typeof interaction.click === "string") {
								if (!aggregatedSelections.has(interaction.click)) {
									aggregatedSelections.set(interaction.click, []);
								}
								aggregatedSelections.get(interaction.click)!.push(numericValue);
							}
						},
					);
				}
			});

			// Apply aggregated selections
			if (stepFrame !== undefined) {
				setStep(stepFrame);
			}

			// Send each geometry selection to server
			aggregatedSelections.forEach((indices, geometryName) => {
				updateSelectionForGeometry(geometryName, indices);
			});
		},
		[setStep, updateSelectionForGeometry],
	);

	const onPlotSelected = useCallback(
		(data: PlotlyJS.PlotSelectionEvent) => {
			if (!data?.points) {
				return;
			}
			const aggregatedSelections = new Map<string, number[]>();
			const selectedFrames: number[] = [];

			data.points.forEach((point) => {
				if (
					point.data?.meta?.interactions &&
					Array.isArray(point.data.meta.interactions)
				) {
					// For each interaction dimension, extract the corresponding customdata value
					point.data.meta.interactions.forEach(
						(
							interaction:
								| { click?: string; select?: string; hover?: string }
								| null
								| undefined,
							dimensionIndex: number,
						) => {
							// Skip if interaction is null/undefined
							if (!interaction) return;

							// Skip if interaction has no select action
							if (!interaction.select) return;

							const numericValue = getCustomDataValue(
								point.customdata,
								dimensionIndex,
							);

							// Skip if customdata value is null/undefined
							if (numericValue === undefined) return;

							// Handle "step" separately for frame selection
							if (interaction.select === "step") {
								selectedFrames.push(numericValue);
							} else if (typeof interaction.select === "string") {
								// Aggregate geometry selections
								if (!aggregatedSelections.has(interaction.select)) {
									aggregatedSelections.set(interaction.select, []);
								}
								aggregatedSelections
									.get(interaction.select)!
									.push(numericValue);
							}
						},
					);
				}
			});

			// Apply frame selection - replace with newly selected frames (dedupe to avoid redundant API calls)
			if (selectedFrames.length > 0) {
				updateFrameSelection([...new Set(selectedFrames)]);
			}

			// Send each geometry selection to server
			aggregatedSelections.forEach((indices, geometryName) => {
				updateSelectionForGeometry(geometryName, indices);
			});
		},
		[updateFrameSelection, updateSelectionForGeometry],
	);

	const onPlotDeselect = useCallback(() => {
		// Guard: avoid redundant API calls if already empty
		if (frame_selection && frame_selection.length > 0) {
			updateFrameSelection([]);
		}
	}, [frame_selection, updateFrameSelection]);

	// ===== EFFECT: Initialize/Update Plotly Chart =====
	useEffect(() => {
		if (!plotContainer.current || !plotlyJson) {
			return;
		}

		// ===== BUILD MARKER TRACKS =====
		// Extract ALL interactions (click and select) to create marker tracks
		const newMarkerTracks: MarkerTrackInfo[] = [];
		const baseDataLength = plotlyJson.data.length;

		plotlyJson.data.forEach((trace: any, dataIdx: number) => {
			if (!trace.meta?.interactions) return;

			const xConverted = convertToArray(trace.x);
			const yConverted = convertToArray(trace.y);
			const customdataConverted = trace.customdata
				? convertToArray(trace.customdata)
				: [];

			// Process each dimension in the interactions array
			trace.meta.interactions.forEach(
				(interaction: any, dimensionIdx: number) => {
					if (!interaction) return; // Skip null interactions

					// Handle click: "step" (frame jumping)
					if (interaction.click === "step") {
						newMarkerTracks.push({
							dataTraceIndex: dataIdx,
							markerTraceIndex: baseDataLength + newMarkerTracks.length,
							xData: xConverted,
							yData: yConverted,
							customdata: customdataConverted,
							interactionType: "step",
							interactionKey: "step",
							customdataIndex: dimensionIdx,
						});
					}

					// Handle select: "step" or select: "<geometry>" - all go into one "selection" marker type
					if (
						interaction.select === "step" ||
						(typeof interaction.select === "string" &&
							interaction.select !== "step")
					) {
						newMarkerTracks.push({
							dataTraceIndex: dataIdx,
							markerTraceIndex: baseDataLength + newMarkerTracks.length,
							xData: xConverted,
							yData: yConverted,
							customdata: customdataConverted,
							interactionType: "selection",
							interactionKey: interaction.select || "step",
							customdataIndex: dimensionIdx,
						});
					}
				},
			);
		});

		// ===== BUILD MARKER TRACES =====
		// Create one marker trace per interaction
		const markerTraces = newMarkerTracks.map((track) => {
			const config = MARKER_CONFIGS[track.interactionType];

			// Initialize markers based on current state
			let initialPoints: Array<{ x: number; y: number }> = [];
			const initialText: string[] = [];

			if (track.interactionType === "step") {
				const point = findPointAtFrame(track, currentFrame);
				if (point) {
					initialPoints = [point];
					initialText.push(`<b>Current Frame</b><br>Frame: ${currentFrame}`);
				}
			} else if (track.interactionType === "selection") {
				// Initialize with existing selections from store
				// Handle frame_selection (interactionKey === "step")
				if (
					track.interactionKey === "step" &&
					frame_selection &&
					Array.isArray(frame_selection) &&
					frame_selection.length > 0
				) {
					frame_selection.forEach((frame) => {
						if (frame === null || frame === undefined) return;
						const points = findPointsWithValue(
							track,
							frame,
							track.customdataIndex,
						);
						points.forEach((p) => {
							initialPoints.push(p);
							initialText.push(`<b>Selected Frame</b><br>Frame: ${frame}`);
						});
					});
				} else if (track.interactionKey !== "step") {
					// Handle geometry selections (arbitrary dimensions)
					const geometryName = track.interactionKey;
					const selectedIndices = selections[geometryName];
					if (Array.isArray(selectedIndices) && selectedIndices.length > 0) {
						selectedIndices.forEach((index) => {
							if (index === null || index === undefined) return;
							const points = findPointsWithValue(
								track,
								index,
								track.customdataIndex,
							);
							points.forEach((p) => {
								initialPoints.push(p);
								initialText.push(
									`<b>Selected ${geometryName.charAt(0).toUpperCase() + geometryName.slice(1)}</b><br>ID: ${index}`,
								);
							});
						});
					}
				}
			}

			return {
				x: initialPoints.map((p) => p.x),
				y: initialPoints.map((p) => p.y),
				text: initialText,
				mode: "markers",
				type: "scatter",
				marker: {
					size: config.size,
					color: config.color,
					line: {
						color: config.borderColor,
						width: config.borderWidth,
					},
				},
				showlegend: false,
				hoverinfo: "text",
				hovertemplate: "%{text}<extra></extra>",
				name: `_marker_${track.interactionType}_${track.interactionKey}_${track.customdataIndex}`,
			};
		});

		// ===== PREPARE FINAL PLOT DATA =====
		// Sort marker traces so "step" markers (current frame) always render on top of "selection" markers
		// Plotly renders traces in order, so later traces appear on top
		const selectionMarkers = markerTraces.filter((t) =>
			t.name?.includes("_marker_selection"),
		);
		const stepMarkers = markerTraces.filter((t) =>
			t.name?.includes("_marker_step"),
		);
		const sortedMarkerTraces = [...selectionMarkers, ...stepMarkers];

		// Update markerTracks with new indices after reordering
		const markerTraceIndexMap = new Map<string, number>();
		sortedMarkerTraces.forEach((trace, idx) => {
			markerTraceIndexMap.set(trace.name || "", baseDataLength + idx);
		});

		const updatedMarkerTracks = newMarkerTracks.map((track) => {
			const nameKey = `_marker_${track.interactionType}_${track.interactionKey}_${track.customdataIndex}`;
			return {
				...track,
				markerTraceIndex:
					markerTraceIndexMap.get(nameKey) || track.markerTraceIndex,
			};
		});
		markerTracks.current = updatedMarkerTracks;

		const finalPlotData = [...plotlyJson.data, ...sortedMarkerTraces];

		// Use Plotly.react for efficient updates
		// This only re-renders changed traces and preserves interaction state
		Plotly.react(plotContainer.current, finalPlotData, plotLayout, {
			responsive: true,
			displaylogo: false,
			displayModeBar: true,
		});

		const container = plotContainer.current;

		// Store handler references for cleanup (Plotly extends HTMLElement at runtime)
		onClickRef.current = onPlotClick;
		onSelectedRef.current = onPlotSelected;
		onDeselectRef.current = onPlotDeselect;

		// Attach event listeners using Plotly.js API
		container.on("plotly_click", onPlotClick);
		container.on("plotly_selected", onPlotSelected);
		container.on("plotly_deselect", onPlotDeselect);

		// Cleanup: Remove listeners by calling removeAllListeners with event name
		return () => {
			if (container) {
				container.removeAllListeners("plotly_click");
				container.removeAllListeners("plotly_selected");
				container.removeAllListeners("plotly_deselect");
			}
		};
	}, [
		plotlyJson,
		plotLayout,
		onPlotClick,
		onPlotSelected,
		onPlotDeselect,
		frame_selection,
		selections,
		currentFrame,
	]);

	// ===== EFFECT: Purge Plotly on unmount =====
	// Read plotContainer.current inside the cleanup so we capture the live ref
	// at unmount, not the (possibly null) ref at first mount.
	useEffect(() => {
		return () => {
			const container = plotContainer.current;
			if (container) {
				Plotly.purge(container);
			}
		};
	}, []);

	// ===== EFFECT: Re-layout Plotly on container size changes =====
	// Dockview emits onDidDimensionsChange on splitter drag, popout, maximize,
	// and viewport resize. Plotly's responsive:true handles window resizes but
	// not intra-panel resizes, so we re-trigger its layout pass here.
	useEffect(() => {
		const disposable = props.api.onDidDimensionsChange(() => {
			const container = plotContainer.current;
			if (container) {
				Plotly.Plots.resize(container);
			}
		});
		return () => disposable.dispose();
	}, [props.api]);

	// ===== EFFECT: Update Step Markers on Frame Change =====
	// Updates red markers for click="step" interaction
	useEffect(() => {
		if (!plotContainer.current || markerTracks.current.length === 0) {
			return;
		}

		const stepMarkers = markerTracks.current.filter(
			(t) => t.interactionType === "step",
		);
		if (stepMarkers.length === 0) return;

		const xUpdates: (number[] | [])[] = [];
		const yUpdates: (number[] | [])[] = [];
		const textUpdates: (string[] | [])[] = [];
		const indices: number[] = [];

		stepMarkers.forEach((track) => {
			indices.push(track.markerTraceIndex);

			const point = findPointAtFrame(track, currentFrame);
			if (point) {
				xUpdates.push([point.x]);
				yUpdates.push([point.y]);
				textUpdates.push([`<b>Current Frame</b><br>Frame: ${currentFrame}`]);
			} else {
				xUpdates.push([]);
				yUpdates.push([]);
				textUpdates.push([]);
			}
		});

		// EFFICIENT UPDATE: Only restyle marker positions and text, no full re-render!
		Plotly.restyle(
			plotContainer.current,
			{ x: xUpdates, y: yUpdates, text: textUpdates as any },
			indices,
		);
	}, [currentFrame]);

	// ===== EFFECT: Update Selection Markers =====
	// Combines and updates all selection markers (frame_selection + geometry selections)
	useEffect(() => {
		if (!plotContainer.current || markerTracks.current.length === 0) {
			return;
		}

		const selectionMarkers = markerTracks.current.filter(
			(t) => t.interactionType === "selection",
		);
		if (selectionMarkers.length === 0) return;

		const xUpdates: (number[] | [])[] = [];
		const yUpdates: (number[] | [])[] = [];
		const textUpdates: (string[] | [])[] = [];
		const indices: number[] = [];

		selectionMarkers.forEach((track) => {
			indices.push(track.markerTraceIndex);

			const allPoints: Array<{ x: number; y: number }> = [];
			const tooltips: string[] = [];

			// Handle frame_selection (interactionKey === "step")
			if (
				track.interactionKey === "step" &&
				Array.isArray(frame_selection) &&
				frame_selection.length > 0
			) {
				frame_selection.forEach((frame) => {
					if (frame === null || frame === undefined) return;
					const points = findPointsWithValue(
						track,
						frame,
						track.customdataIndex,
					);
					points.forEach((p) => {
						allPoints.push(p);
						tooltips.push(`<b>Selected Frame</b><br>Frame: ${frame}`);
					});
				});
			} else if (track.interactionKey !== "step") {
				// Handle geometry selections (arbitrary dimensions)
				const geometryName = track.interactionKey;
				const selectedIndices = selections[geometryName];

				if (Array.isArray(selectedIndices) && selectedIndices.length > 0) {
					selectedIndices.forEach((index) => {
						if (index === null || index === undefined) return;
						const points = findPointsWithValue(
							track,
							index,
							track.customdataIndex,
						);
						points.forEach((p) => {
							allPoints.push(p);
							tooltips.push(
								`<b>Selected ${geometryName.charAt(0).toUpperCase() + geometryName.slice(1)}</b><br>ID: ${index}`,
							);
						});
					});
				}
			}

			xUpdates.push(allPoints.map((p) => p.x));
			yUpdates.push(allPoints.map((p) => p.y));
			textUpdates.push(tooltips);
		});

		Plotly.restyle(
			plotContainer.current,
			{ x: xUpdates, y: yUpdates, text: textUpdates as any },
			indices,
		);
	}, [frame_selection, selections]);

	// ===== RENDER: Loading State =====
	if (isLoading) {
		return (
			<Box
				data-testid={`plot-view-${figureKey}`}
				sx={{
					width: "100%",
					height: "100%",
					display: "flex",
					justifyContent: "center",
					alignItems: "center",
					backgroundColor: "background.paper",
				}}
			>
				<CircularProgress />
			</Box>
		);
	}

	// ===== RENDER: Error State =====
	if (isError) {
		return (
			<Box
				data-testid={`plot-view-${figureKey}`}
				sx={{
					width: "100%",
					height: "100%",
					display: "flex",
					flexDirection: "column",
					justifyContent: "center",
					alignItems: "center",
					backgroundColor: "background.paper",
				}}
			>
				<ErrorIcon sx={{ fontSize: 40, mb: 1, color: "error.main" }} />
				<Typography variant="body1" color="error">
					Error loading figure.
				</Typography>
				<Typography variant="caption" color="error">
					{error?.message}
				</Typography>
			</Box>
		);
	}

	// ===== RENDER: Invalid Figure Type =====
	if (!figureResponse?.figure) {
		return (
			<Box
				data-testid={`plot-view-${figureKey}`}
				sx={{
					width: "100%",
					height: "100%",
					p: 1,
					backgroundColor: "background.paper",
				}}
			>
				<Typography color="error">No figure data available.</Typography>
			</Box>
		);
	}

	// ===== RENDER: Plotly Figure =====
	if (figureResponse.figure.type === "plotly") {
		return (
			<Box
				data-testid={`plot-view-${figureKey}`}
				ref={plotContainer}
				sx={{
					width: "100%",
					height: "100%",
					overflow: "auto",
					position: "relative",
					backgroundColor: "background.paper",
				}}
			/>
		);
	}

	// ===== RENDER: Other Figure Types =====
	return (
		<Box
			data-testid={`plot-view-${figureKey}`}
			sx={{
				width: "100%",
				height: "100%",
				p: 1,
				overflow: "auto",
				backgroundColor: "background.paper",
			}}
		>
			<pre style={{ margin: 0 }}>
				{JSON.stringify(figureResponse.figure, null, 2)}
			</pre>
		</Box>
	);
}
