import React, { useRef, useEffect, useCallback, useMemo, useState } from "react";
import { Rnd, Position } from "react-rnd";
import Plotly from "plotly.js-dist-min";
import type * as PlotlyJS from "plotly.js-dist-min";
import { decodeTypedArraySpec } from "plotly.js/src/lib/array.js";

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
import { useWindowManagerStore } from "../stores/windowManagerStore";
import { useFigure, useFigureList } from "../hooks/useFigures";
import { useAppStore } from "../store";
import { useStepControl } from "../hooks/useStepControl";

// MUI Imports
import {
  Box,
  Paper,
  Typography,
  IconButton,
  CircularProgress,
  TextField,
  Autocomplete,
  Tooltip,
} from "@mui/material";
import {
  Close as CloseIcon,
  ErrorOutline as ErrorIcon,
} from "@mui/icons-material";

interface FigureWindowProps {
  windowId: string;
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
 * Window template component that eliminates UI boilerplate duplication.
 * Provides consistent window structure across all states (loading, error, content).
 */
interface WindowTemplateProps {
  windowInstance: any;
  allFiguresResponse: any;
  handleDragStart: () => void;
  handleDragStop: (event: any, d: Position) => void;
  handleResizeStop: (event: any, direction: any, ref: HTMLElement, delta: any, position: Position) => void;
  handleMouseDown: () => void;
  handleClose: () => void;
  handleAutocompleteChange: (event: React.SyntheticEvent, newValue: string | null) => void;
  markerControls?: React.ReactNode; // Optional marker visibility controls
  children: React.ReactNode;
}

const WindowTemplate = React.memo(
  ({
    windowInstance,
    allFiguresResponse,
    handleDragStart,
    handleDragStop,
    handleResizeStop,
    handleMouseDown,
    handleClose,
    handleAutocompleteChange,
    markerControls,
    children,
  }: WindowTemplateProps) => (
    <Rnd
      size={{ width: windowInstance.width, height: windowInstance.height }}
      position={{ x: windowInstance.x, y: windowInstance.y }}
      minWidth={350}
      minHeight={300}
      style={{ zIndex: windowInstance.zIndex }}
      dragHandleClassName="drag-handle"
      onDragStart={handleDragStart}
      onDragStop={handleDragStop}
      onResizeStop={handleResizeStop}
      bounds=".drag-boundary-container"
    >
      <Paper
        elevation={4}
        sx={{
          width: "100%",
          height: "100%",
          display: "flex",
          flexDirection: "column",
          overflow: "hidden",
        }}
      >
        <Box
          className="drag-handle"
          onMouseDown={handleMouseDown}
          sx={{
            p: 1,
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
            borderBottom: 1,
            borderColor: "divider",
            cursor: "move",
            gap: 1,
          }}
        >
          <Autocomplete
            options={allFiguresResponse?.figures || []}
            value={windowInstance.figureKey}
            onChange={handleAutocompleteChange}
            disableClearable
            size="small"
            onMouseDown={(event) => event.stopPropagation()}
            disablePortal
            sx={{ flexGrow: 1, minWidth: 150 }}
            renderInput={(params) => (
              <TextField
                {...params}
                label="Figure"
                variant="outlined"
                size="small"
              />
            )}
          />
          {markerControls}
          <IconButton onClick={handleClose} size="small">
            <CloseIcon />
          </IconButton>
        </Box>

        {children}
      </Paper>
    </Rnd>
  ),
);

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
        if (flatArray && typeof flatArray === "object" && flatArray.length !== undefined) {
          return Array.from(flatArray);
        }
        return [];
      }

      // If data has a shape parameter, reshape the array to match dimensions
      if (data.shape) {
        const shapeStr = data.shape;
        // Parse shape like "67, 2" to get all dimensions
        const dims = shapeStr.split(",").map((s: string) => parseInt(s.trim()));

        if (dims.length >= 1) {
          // Recursively reshape for arbitrary dimensions
          const result = reshapeFlat(flatArray, dims, 0);
          return result;
        }
      }

      return flatArray;
    } catch (e) {
      return [];
    }
  }

  return [];
}

/**
 * Recursively reshape a flat array into an n-dimensional array.
 * @param flatArray - The flat array to reshape
 * @param dims - Array of dimensions [d0, d1, d2, ...]
 * @param dimIndex - Current dimension being processed
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
    const subArray = reshapeFlat(flatArray.slice(start, end), dims, dimIndex + 1);
    result.push(subArray);
  }

  return result;
}

/**
 * [FIXED] Extracts a single numeric value from a customdata entry at a specific dimension.
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
      const numKeys = Object.keys(element).filter((k) => !isNaN(Number(k)));
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
    const numKeys = Object.keys(customDataPoint).filter((k) => !isNaN(Number(k)));
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
  return isNaN(numericValue) ? undefined : numericValue;
}

/**
 * Find data points matching a specific value in a customdata dimension.
 * Returns array of {x, y} positions for all matching points.
 *
 * [REFACTORED] Uses getCustomDataValue helper.
 */
function findPointsWithValue(
  track: MarkerTrackInfo,
  value: number,
  dimensionIndex: number = 0,
): Array<{ x: number; y: number }> {
  const points: Array<{ x: number; y: number }> = [];

  // Case 1: customdata exists - use it to find matching values
  if (track.customdata.length > 0) {
    for (let i = 0; i < track.customdata.length; i++) {
      const customDataPoint = track.customdata[i];

      // Use the new helper function
      const pointValue = getCustomDataValue(customDataPoint, dimensionIndex);

      // Use loose equality to handle string vs number comparison
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
 * Kept for backwards compatibility with existing step marker logic.
 */
function findPointAtFrame(
  track: MarkerTrackInfo,
  frame: number,
): { x: number; y: number } | null {
  const points = findPointsWithValue(track, frame, track.customdataIndex);
  return points.length > 0 ? points[0] : null;
}

function FigureWindow({ windowId }: FigureWindowProps) {
  // ===== LOCAL STATE: Marker Visibility Toggles =====
  const [markerVisibility, setMarkerVisibility] = useState({
    step: true, // Red markers for current frame
    selection: true, // Blue markers for all selections (frame_selection + geometry selections)
  });

  // ===== DOM REFERENCES =====
  // Plotly extends the HTML element at runtime with event methods
  const plotContainer = useRef<PlotlyJS.PlotlyHTMLElement>(null);

  // ===== MARKER TRACKING =====
  // Stores info about marker traces for efficient frame updates
  const markerTracks = useRef<MarkerTrackInfo[]>([]);

  // ===== EVENT LISTENER REFS =====
  // Keep references to handlers for proper cleanup
  const onClickRef = useRef<((data: PlotlyJS.PlotMouseEvent) => void) | null>(null);
  const onSelectedRef = useRef<((data: PlotlyJS.PlotSelectionEvent) => void) | null>(null);
  const onDeselectRef = useRef<(() => void) | null>(null);

  // ===== STORE SELECTORS =====
  const windowInstance = useWindowManagerStore(
    (state) => state.openWindows[windowId],
  );
  const { updateWindowState, closeWindow, bringToFront, changeFigureInWindow } =
    useWindowManagerStore();

  // ===== DATA FETCHING =====
  const {
    data: figureResponse,
    isLoading,
    isError,
    error,
  } = useFigure(windowInstance?.figureKey, { enabled: !!windowInstance });
  const { data: allFiguresResponse } = useFigureList();

  // ===== STORE ACTIONS & STATE =====
  // Use individual selectors to prevent unnecessary re-renders
  const updateSelectionForGeometry = useAppStore((state) => state.updateSelectionForGeometry);
  const setFrameSelection = useAppStore((state) => state.setFrameSelection);
  const currentFrame = useAppStore((state) => state.currentFrame);
  const frame_selection = useAppStore((state) => state.frame_selection);
  const selections = useAppStore((state) => state.selections);
  const { setStep } = useStepControl();

  if (!windowInstance) {
    return null;
  }

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
          figureKey: windowInstance?.figureKey,
          error: e instanceof Error ? e.message : String(e),
        });
        return null;
      }
    }
    return null;
  }, [figureResponse?.figure, windowInstance?.figureKey]);

  // ===== MEMOIZED: Layout with autosize =====
  const plotLayout = useMemo(
    () => ({
      ...plotlyJson?.layout,
      autosize: true,
    }),
    [plotlyJson?.layout],
  );

  // ===== STABLE CALLBACKS: Event Handlers =====

  const onPlotClick = useCallback(
    (data: PlotlyJS.PlotMouseEvent) => {
      const aggregatedSelections = new Map<string, number[]>();
      let stepFrame: number | undefined;

      data.points.forEach((point) => {
        if (point.data?.meta?.interactions && Array.isArray(point.data.meta.interactions)) {
          // For each interaction dimension, extract the corresponding customdata value
          point.data.meta.interactions.forEach(
            (interaction: { click?: string; select?: string; hover?: string } | null | undefined, dimensionIndex: number) => {
              // Skip if interaction is null/undefined
              if (!interaction) return;

              // Skip if interaction has no click action
              if (!interaction.click) return;

              // [FIXED] Extract value using the robust helper function
              const numericValue = getCustomDataValue(point.customdata, dimensionIndex);

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
      let selectedFrames: number[] = [];

      data.points.forEach((point) => {
        if (point.data?.meta?.interactions && Array.isArray(point.data.meta.interactions)) {
          // For each interaction dimension, extract the corresponding customdata value
          point.data.meta.interactions.forEach(
            (interaction: { click?: string; select?: string; hover?: string } | null | undefined, dimensionIndex: number) => {
              // Skip if interaction is null/undefined
              if (!interaction) return;

              // Skip if interaction has no select action
              if (!interaction.select) return;

              // [FIXED] Extract value using the robust helper function
              const numericValue = getCustomDataValue(point.customdata, dimensionIndex);

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
                aggregatedSelections.get(interaction.select)!.push(numericValue);
              }
            },
          );
        }
      });

      // Apply frame selection - replace with newly selected frames
      if (selectedFrames.length > 0) {
        setFrameSelection(selectedFrames);
      }

      // Send each geometry selection to server
      aggregatedSelections.forEach((indices, geometryName) => {
        updateSelectionForGeometry(geometryName, indices);
      });
    },
    [setFrameSelection, updateSelectionForGeometry],
  );

  const onPlotDeselect = useCallback(() => {
    setFrameSelection([]);
  }, [setFrameSelection]);

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
      trace.meta.interactions.forEach((interaction: any, dimensionIdx: number) => {
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
        if (interaction.select === "step" || (typeof interaction.select === "string" && interaction.select !== "step")) {
          newMarkerTracks.push({
            dataTraceIndex: dataIdx,
            markerTraceIndex: baseDataLength + newMarkerTracks.length,
            xData: xConverted,
            yData: yConverted,
            customdata: customdataConverted,
            interactionType: "selection",
            interactionKey: interaction.select || "step", // Use the select value as key
            customdataIndex: dimensionIdx,
          });
        }
      });
    });

    // ===== BUILD MARKER TRACES =====
    // Create one marker trace per interaction
    const markerTraces = newMarkerTracks.map((track) => {
      const config = MARKER_CONFIGS[track.interactionType];

      // Initialize markers based on current state
      let initialPoints: Array<{ x: number; y: number }> = [];
      let initialText: string[] = [];

      if (track.interactionType === "step") {
        const point = findPointAtFrame(track, currentFrame);
        if (point) {
          initialPoints = [point];
          initialText = [`<b>Current Frame</b><br>Frame: ${currentFrame}`];
        }
      } else if (track.interactionType === "selection") {
        // Initialize with existing selections from store
        // Handle frame_selection (interactionKey === "step")
        if (track.interactionKey === "step" && frame_selection && Array.isArray(frame_selection) && frame_selection.length > 0) {
          frame_selection.forEach((frame) => {
            // Skip null/undefined frames
            if (frame === null || frame === undefined) return;
            const points = findPointsWithValue(track, frame, track.customdataIndex);
            points.forEach((p) => {
              initialPoints.push(p);
              initialText.push(`<b>Selected Frame</b><br>Frame: ${frame}`);
            });
          });
        }
        // Handle geometry selections (arbitrary dimensions)
        else if (track.interactionKey !== "step") {
          const geometryName = track.interactionKey;
          const selectedIndices = selections[geometryName];
          if (Array.isArray(selectedIndices) && selectedIndices.length > 0) {
            selectedIndices.forEach((index) => {
              // Skip null/undefined indices
              if (index === null || index === undefined) return;
              const points = findPointsWithValue(track, index, track.customdataIndex);
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
    const selectionMarkers = markerTraces.filter((t) => t.name?.includes("_marker_selection"));
    const stepMarkers = markerTraces.filter((t) => t.name?.includes("_marker_step"));
    const sortedMarkerTraces = [...selectionMarkers, ...stepMarkers];

    // Update markerTracks with new indices after reordering
    const markerTraceIndexMap = new Map<string, number>();
    sortedMarkerTraces.forEach((trace, idx) => {
      markerTraceIndexMap.set(trace.name || "", baseDataLength + idx);
    });

    // [FIXED] Use the correct, full name (including customdataIndex) to find the new trace index
    const updatedMarkerTracks = newMarkerTracks.map((track) => {
      const nameKey = `_marker_${track.interactionType}_${track.interactionKey}_${track.customdataIndex}`;
      return {
        ...track,
        markerTraceIndex: markerTraceIndexMap.get(nameKey) || track.markerTraceIndex,
      };
    });
    markerTracks.current = updatedMarkerTracks;

    const finalPlotData = [...plotlyJson.data, ...sortedMarkerTraces];

    // Use Plotly.react for efficient updates
    // This only re-renders changed traces and preserves interaction state
    Plotly.react(plotContainer.current, finalPlotData, plotLayout, {
      responsive: true,
      displayModeBar: true,
    });

    const container = plotContainer.current;

    // Store handler references for cleanup (Plotly extends HTMLElement at runtime)
    onClickRef.current = onPlotClick;
    onSelectedRef.current = onPlotSelected;
    onDeselectRef.current = onPlotDeselect;

    // Attach event listeners using Plotly.js API
    // Plotly.js provides .on() method on the graph div element
    container.on("plotly_click", onPlotClick);
    container.on("plotly_selected", onPlotSelected);
    container.on("plotly_deselect", onPlotDeselect);

    // Cleanup: Remove listeners by calling removeAllListeners with event name
    // Plotly.js API: removeAllListeners(eventName) removes all listeners for that specific event
    return () => {
      if (container) {
        container.removeAllListeners("plotly_click");
        container.removeAllListeners("plotly_selected");
        container.removeAllListeners("plotly_deselect");
      }
    };
  }, [plotlyJson, plotLayout, onPlotClick, onPlotSelected, frame_selection, selections, currentFrame]);

  // ===== EFFECT: Update Step Markers on Frame Change =====
  // Updates red markers for click="step" interaction
  useEffect(() => {
    if (!plotContainer.current || markerTracks.current.length === 0) {
      return;
    }

    const stepMarkers = markerTracks.current.filter((t) => t.interactionType === "step");
    if (stepMarkers.length === 0) return;

    const xUpdates: (number[] | [])[] = [];
    const yUpdates: (number[] | [])[] = [];
    const textUpdates: (string[] | [])[] = [];
    const indices: number[] = [];

    stepMarkers.forEach((track) => {
      indices.push(track.markerTraceIndex);

      // If visibility is off - hide markers
      if (!markerVisibility.step) {
        xUpdates.push([]);
        yUpdates.push([]);
        textUpdates.push([]);
        return;
      }

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
    Plotly.restyle(plotContainer.current, { x: xUpdates, y: yUpdates, text: textUpdates as any }, indices);
  }, [currentFrame, markerVisibility.step]);

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

      // If visibility is off - hide markers
      if (!markerVisibility.selection) {
        xUpdates.push([]);
        yUpdates.push([]);
        textUpdates.push([]);
        return;
      }

      const allPoints: Array<{ x: number; y: number }> = [];
      const tooltips: string[] = [];

      // Handle frame_selection (interactionKey === "step")
      if (track.interactionKey === "step" && Array.isArray(frame_selection) && frame_selection.length > 0) {
        frame_selection.forEach((frame) => {
          // Skip null/undefined frames
          if (frame === null || frame === undefined) return;
          const points = findPointsWithValue(track, frame, track.customdataIndex);
          points.forEach((p) => {
            allPoints.push(p);
            tooltips.push(`<b>Selected Frame</b><br>Frame: ${frame}`);
          });
        });
      }
      // Handle geometry selections (arbitrary dimensions)
      else if (track.interactionKey !== "step") {
        const geometryName = track.interactionKey;
        const selectedIndices = selections[geometryName];

        if (Array.isArray(selectedIndices) && selectedIndices.length > 0) {
          selectedIndices.forEach((index) => {
            // Skip null/undefined indices
            if (index === null || index === undefined) return;
            const points = findPointsWithValue(track, index, track.customdataIndex);
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

    Plotly.restyle(plotContainer.current, { x: xUpdates, y: yUpdates, text: textUpdates as any }, indices);
  }, [frame_selection, selections, markerVisibility.selection]);

  // ===== STABLE CALLBACKS: Window Management =====

  const handleAutocompleteChange = useCallback(
    (_event: React.SyntheticEvent, newValue: string | null) => {
      if (newValue) {
        changeFigureInWindow(windowId, newValue);
      }
    },
    [windowId, changeFigureInWindow],
  );

  const handleDragStart = useCallback(() => {
    bringToFront(windowId);
  }, [windowId, bringToFront]);

  const handleDragStop = useCallback(
    (_event: any, d: Position) => {
      updateWindowState(windowId, { x: d.x, y: d.y });
    },
    [windowId, updateWindowState],
  );

  const handleResizeStop = useCallback(
    (_event: any, _direction: any, ref: HTMLElement, _delta: any, position: Position) => {
      const newWidth = ref.offsetWidth;
      const newHeight = ref.offsetHeight;

      // Update window state
      updateWindowState(windowId, {
        width: newWidth,
        height: newHeight,
        ...position,
      });

      // Trigger Plotly resize via relayout (only updates layout, preserves data and markers)
      if (plotContainer.current) {
        Plotly.relayout(plotContainer.current, plotLayout);
      }
    },
    [windowId, updateWindowState, plotLayout],
  );

  const handleMouseDown = useCallback(() => {
    bringToFront(windowId);
  }, [windowId, bringToFront]);

  const handleClose = useCallback(() => {
    closeWindow(windowId);
  }, [windowId, closeWindow]);

  // ===== MARKER VISIBILITY CONTROLS =====
  // Render toggle buttons for marker visibility (only for Plotly figures)
  const markerControls = useMemo(() => {
    if (figureResponse?.figure?.type !== "plotly") return null;

    return (
      <Box sx={{ display: "flex", gap: 0.5, alignItems: "center" }}>
        <Typography variant="caption" sx={{ fontSize: "0.7rem", mr: 0.5 }}>
          Markers:
        </Typography>

        {/* Frame Marker Toggle */}
        <Tooltip
          title={markerVisibility.step ? "Hide frame marker (red)" : "Show frame marker (red)"}
          placement="bottom"
        >
          <IconButton
            size="small"
            onClick={() =>
              setMarkerVisibility((prev) => ({ ...prev, step: !prev.step }))
            }
            sx={{
              p: 0.25,
              bgcolor: markerVisibility.step ? "rgba(255, 0, 0, 0.1)" : "transparent",
              borderRadius: 1,
              "&:hover": { bgcolor: markerVisibility.step ? "rgba(255, 0, 0, 0.2)" : "rgba(0, 0, 0, 0.05)" },
            }}
          >
            <Box
              sx={{
                width: 8,
                height: 8,
                borderRadius: "50%",
                backgroundColor: "red",
              }}
            />
          </IconButton>
        </Tooltip>

        {/* Selection Markers Toggle (combines frame selection + geometry selections) */}
        <Tooltip
          title={markerVisibility.selection ? "Hide selection markers (blue)" : "Show selection markers (blue)"}
          placement="bottom"
        >
          <IconButton
            size="small"
            onClick={() =>
              setMarkerVisibility((prev) => ({ ...prev, selection: !prev.selection }))
            }
            sx={{
              p: 0.25,
              bgcolor: markerVisibility.selection ? "rgba(70, 130, 180, 0.1)" : "transparent",
              borderRadius: 1,
              "&:hover": { bgcolor: markerVisibility.selection ? "rgba(70, 130, 180, 0.2)" : "rgba(0, 0, 0, 0.05)" },
            }}
          >
            <Box
              sx={{
                width: 8,
                height: 8,
                borderRadius: "50%",
                backgroundColor: "rgb(70, 130, 180)",
              }}
            />
          </IconButton>
        </Tooltip>
      </Box>
    );
  }, [figureResponse?.figure?.type, markerVisibility]);

  // ===== RENDER: Loading State =====
  if (isLoading) {
    return (
      <WindowTemplate
        windowInstance={windowInstance}
        allFiguresResponse={allFiguresResponse}
        handleDragStart={handleDragStart}
        handleDragStop={handleDragStop}
        handleResizeStop={handleResizeStop}
        handleMouseDown={handleMouseDown}
        handleClose={handleClose}
        handleAutocompleteChange={handleAutocompleteChange}
        markerControls={markerControls}
      >
        <Box
          sx={{
            flexGrow: 1,
            p: 1,
            overflow: "auto",
            position: "relative",
            backgroundColor: "background.paper",
            display: "flex",
            justifyContent: "center",
            alignItems: "center",
          }}
        >
          <CircularProgress />
        </Box>
      </WindowTemplate>
    );
  }

  // ===== RENDER: Error State =====
  if (isError) {
    return (
      <WindowTemplate
        windowInstance={windowInstance}
        allFiguresResponse={allFiguresResponse}
        handleDragStart={handleDragStart}
        handleDragStop={handleDragStop}
        handleResizeStop={handleResizeStop}
        handleMouseDown={handleMouseDown}
        handleClose={handleClose}
        handleAutocompleteChange={handleAutocompleteChange}
        markerControls={markerControls}
      >
        <Box
          sx={{
            flexGrow: 1,
            p: 1,
            overflow: "auto",
            position: "relative",
            backgroundColor: "background.paper",
            display: "flex",
            flexDirection: "column",
            justifyContent: "center",
            alignItems: "center",
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
      </WindowTemplate>
    );
  }

  // ===== RENDER: Invalid Figure Type =====
  if (!figureResponse?.figure) {
    return (
      <WindowTemplate
        windowInstance={windowInstance}
        allFiguresResponse={allFiguresResponse}
        handleDragStart={handleDragStart}
        handleDragStop={handleDragStop}
        handleResizeStop={handleResizeStop}
        handleMouseDown={handleMouseDown}
        handleClose={handleClose}
        handleAutocompleteChange={handleAutocompleteChange}
        markerControls={markerControls}
      >
        <Box
          sx={{
            flexGrow: 1,
            p: 1,
            overflow: "auto",
            position: "relative",
            backgroundColor: "background.paper",
          }}
        >
          <Typography color="error">No figure data available.</Typography>
        </Box>
      </WindowTemplate>
    );
  }

  // ===== RENDER: Plotly Figure =====
  if (figureResponse.figure.type === "plotly" && !isLoading) {
    return (
      <WindowTemplate
        windowInstance={windowInstance}
        allFiguresResponse={allFiguresResponse}
        handleDragStart={handleDragStart}
        handleDragStop={handleDragStop}
        handleResizeStop={handleResizeStop}
        handleMouseDown={handleMouseDown}
        handleClose={handleClose}
        handleAutocompleteChange={handleAutocompleteChange}
        markerControls={markerControls}
      >
        {/* Plotly Chart Container - Direct DOM manipulation via Plotly.react */}
        <Box
          ref={plotContainer}
          sx={{
            flexGrow: 1,
            width: "100%",
            height: "100%",
            overflow: "auto",
            position: "relative",
            backgroundColor: "background.paper",
          }}
        />
      </WindowTemplate>
    );
  }

  // ===== RENDER: Other Figure Types =====
  return (
    <WindowTemplate
      windowInstance={windowInstance}
      allFiguresResponse={allFiguresResponse}
      handleDragStart={handleDragStart}
      handleDragStop={handleDragStop}
      handleResizeStop={handleResizeStop}
      handleMouseDown={handleMouseDown}
      handleClose={handleClose}
      handleAutocompleteChange={handleAutocompleteChange}
      markerControls={markerControls}
    >
      <Box
        sx={{
          flexGrow: 1,
          p: 1,
          overflow: "auto",
          position: "relative",
          backgroundColor: "background.paper",
        }}
      >
        <pre style={{ margin: 0 }}>
          {JSON.stringify(figureResponse.figure, null, 2)}
        </pre>
      </Box>
    </WindowTemplate>
  );
}

export default React.memo(FigureWindow);
