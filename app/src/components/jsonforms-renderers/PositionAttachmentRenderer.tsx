/**
 * PositionAttachmentRenderer - Custom JSONForms control for position fields
 * that support both direct coordinates and CurveAttachment.
 *
 * Renders when schema has x-custom-type="position-attachment"
 *
 * Supports two modes:
 * - Coordinates: Direct [x, y, z] tuple input
 * - Curve Attachment: Reference to a Curve geometry with progress slider
 */

import { useState, useEffect, useCallback } from "react";
import { withJsonFormsControlProps } from "@jsonforms/react";
import {
	rankWith,
	schemaMatches,
	and,
	uiTypeIs,
	ControlProps,
} from "@jsonforms/core";
import {
	Box,
	TextField,
	Typography,
	ToggleButton,
	ToggleButtonGroup,
	Slider,
	FormControl,
	InputLabel,
	Select,
	MenuItem,
	Paper,
	SelectChangeEvent,
} from "@mui/material";
import LocationOnIcon from "@mui/icons-material/LocationOn";
import TimelineIcon from "@mui/icons-material/Timeline";
import { isCurveAttachment, CurveAttachment } from "../../utils/cameraUtils";

type PositionMode = "coordinates" | "curve";

/**
 * Detect if value is a coordinate tuple [x, y, z]
 */
function isCoordinateTuple(
	value: unknown,
): value is [number, number, number] | number[] {
	return (
		Array.isArray(value) &&
		value.length === 3 &&
		value.every((v) => typeof v === "number")
	);
}

/**
 * Get initial mode based on current value
 */
function getInitialMode(value: unknown): PositionMode {
	if (isCurveAttachment(value)) {
		return "curve";
	}
	return "coordinates";
}

/**
 * PositionAttachmentRenderer component
 */
const PositionAttachmentRenderer = ({
	data,
	handleChange,
	path,
	label,
	schema,
	required,
	errors,
	visible,
}: ControlProps) => {
	// Get curve options from injected enum (populated by injectDynamicEnums)
	const curveOptions = (schema.enum as string[]) || [];

	// State for mode toggle
	const [mode, setMode] = useState<PositionMode>(() => getInitialMode(data));

	// Local state for coordinate inputs (to allow typing without immediate validation)
	const [coords, setCoords] = useState<[string, string, string]>(() => {
		if (isCoordinateTuple(data)) {
			return [String(data[0]), String(data[1]), String(data[2])];
		}
		return ["0", "0", "0"];
	});

	// Local state for curve attachment
	const [curveKey, setCurveKey] = useState<string>(() => {
		if (isCurveAttachment(data)) {
			return data.geometry_key;
		}
		return curveOptions[0] || "";
	});

	const [progress, setProgress] = useState<number>(() => {
		if (isCurveAttachment(data)) {
			return data.progress * 100; // Store as percentage for slider
		}
		return 0;
	});

	// Sync local state when external data changes
	useEffect(() => {
		const newMode = getInitialMode(data);
		setMode(newMode);

		if (isCoordinateTuple(data)) {
			setCoords([String(data[0]), String(data[1]), String(data[2])]);
		} else if (isCurveAttachment(data)) {
			setCurveKey(data.geometry_key);
			setProgress(data.progress * 100);
		}
	}, [data]);

	// Update curveKey when curveOptions becomes available (enum injection is async)
	useEffect(() => {
		if (!curveKey && curveOptions.length > 0) {
			setCurveKey(curveOptions[0]);
		}
	}, [curveOptions, curveKey]);

	// Handle mode change
	const handleModeChange = useCallback(
		(_: React.MouseEvent<HTMLElement>, newMode: PositionMode | null) => {
			if (newMode === null) return; // Prevent deselection

			setMode(newMode);

			if (newMode === "coordinates") {
				// Switch to coordinates - use current coords or defaults
				const x = parseFloat(coords[0]) || 0;
				const y = parseFloat(coords[1]) || 0;
				const z = parseFloat(coords[2]) || 0;
				handleChange(path, [x, y, z]);
			} else {
				// Switch to curve attachment
				const attachment: CurveAttachment = {
					type: "curve_attachment",
					geometry_key: curveKey || curveOptions[0] || "",
					progress: progress / 100,
				};
				handleChange(path, attachment);
			}
		},
		[coords, curveKey, progress, curveOptions, handleChange, path],
	);

	// Handle coordinate input changes
	const handleCoordChange = useCallback(
		(index: number, value: string) => {
			const newCoords = [...coords] as [string, string, string];
			newCoords[index] = value;
			setCoords(newCoords);

			// Parse and update if valid
			const x = parseFloat(newCoords[0]);
			const y = parseFloat(newCoords[1]);
			const z = parseFloat(newCoords[2]);

			if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
				handleChange(path, [x, y, z]);
			}
		},
		[coords, handleChange, path],
	);

	// Handle curve selection change
	const handleCurveChange = useCallback(
		(event: SelectChangeEvent<string>) => {
			const newKey = event.target.value;
			setCurveKey(newKey);

			const attachment: CurveAttachment = {
				type: "curve_attachment",
				geometry_key: newKey,
				progress: progress / 100,
			};
			handleChange(path, attachment);
		},
		[progress, handleChange, path],
	);

	// Handle progress slider change
	const handleProgressChange = useCallback(
		(_: Event, newValue: number | number[]) => {
			const newProgress = newValue as number;
			setProgress(newProgress);

			const attachment: CurveAttachment = {
				type: "curve_attachment",
				geometry_key: curveKey,
				progress: newProgress / 100,
			};
			handleChange(path, attachment);
		},
		[curveKey, handleChange, path],
	);

	if (!visible) {
		return null;
	}

	return (
		<Box sx={{ marginBottom: 2 }}>
			<Typography
				variant="subtitle2"
				sx={{ marginBottom: 1, color: "text.secondary" }}
			>
				{label}
				{required && " *"}
			</Typography>

			{/* Mode Toggle */}
			<ToggleButtonGroup
				value={mode}
				exclusive
				onChange={handleModeChange}
				size="small"
				sx={{ marginBottom: 1.5 }}
			>
				<ToggleButton value="coordinates">
					<LocationOnIcon sx={{ mr: 0.5 }} fontSize="small" />
					Coordinates
				</ToggleButton>
				<ToggleButton value="curve" disabled={curveOptions.length === 0}>
					<TimelineIcon sx={{ mr: 0.5 }} fontSize="small" />
					Curve
				</ToggleButton>
			</ToggleButtonGroup>

			<Paper elevation={1} sx={{ padding: 2 }}>
				{mode === "coordinates" ? (
					/* Coordinates Mode */
					<Box sx={{ display: "flex", gap: 1 }}>
						<TextField
							label="X"
							type="number"
							size="small"
							value={coords[0]}
							onChange={(e) => handleCoordChange(0, e.target.value)}
							error={!!errors}
							inputProps={{ step: "any" }}
							sx={{ flex: 1 }}
						/>
						<TextField
							label="Y"
							type="number"
							size="small"
							value={coords[1]}
							onChange={(e) => handleCoordChange(1, e.target.value)}
							error={!!errors}
							inputProps={{ step: "any" }}
							sx={{ flex: 1 }}
						/>
						<TextField
							label="Z"
							type="number"
							size="small"
							value={coords[2]}
							onChange={(e) => handleCoordChange(2, e.target.value)}
							error={!!errors}
							inputProps={{ step: "any" }}
							sx={{ flex: 1 }}
						/>
					</Box>
				) : (
					/* Curve Attachment Mode */
					<Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
						<FormControl fullWidth size="small">
							<InputLabel>Curve Geometry</InputLabel>
							<Select
								value={curveKey}
								label="Curve Geometry"
								onChange={handleCurveChange}
							>
								{curveOptions.length === 0 ? (
									<MenuItem disabled value="">
										No Curve geometries available
									</MenuItem>
								) : (
									curveOptions.map((key) => (
										<MenuItem key={key} value={key}>
											{key}
										</MenuItem>
									))
								)}
							</Select>
						</FormControl>

						<Box>
							<Typography variant="body2" color="text.secondary" gutterBottom>
								Progress: {progress.toFixed(0)}%
							</Typography>
							<Slider
								value={progress}
								onChange={handleProgressChange}
								min={0}
								max={100}
								step={1}
								valueLabelDisplay="auto"
								valueLabelFormat={(v) => `${v}%`}
							/>
						</Box>
					</Box>
				)}
			</Paper>

			{errors && (
				<Typography variant="caption" color="error" sx={{ mt: 0.5 }}>
					{errors}
				</Typography>
			)}
		</Box>
	);
};

/**
 * Tester function for PositionAttachmentRenderer
 * Matches fields with x-custom-type="position-attachment"
 * Priority 10 to override default renderers
 */
export const positionAttachmentTester = rankWith(
	10,
	and(
		schemaMatches(
			(schema) => (schema as any)["x-custom-type"] === "position-attachment",
		),
		uiTypeIs("Control"),
	),
);

export default withJsonFormsControlProps(PositionAttachmentRenderer);
