/**
 * LightPositionRenderer - Custom JSONForms renderer for LightPosition objects.
 *
 * Renders a toggle between "Fixed" (world coordinates) and "Camera" (offset
 * from camera) with three coordinate inputs.
 *
 * Matches when the resolved schema has a "camera_attached" boolean property
 * alongside "x", "y", "z" number properties (i.e. the LightPosition model).
 */

import {
	type ControlProps,
	and,
	rankWith,
	schemaMatches,
	uiTypeIs,
} from "@jsonforms/core";
import { withJsonFormsControlProps } from "@jsonforms/react";
import CameraAltIcon from "@mui/icons-material/CameraAlt";
import LocationOnIcon from "@mui/icons-material/LocationOn";
import {
	Box,
	Paper,
	TextField,
	ToggleButton,
	ToggleButtonGroup,
	Typography,
} from "@mui/material";
import { useCallback, useEffect, useState } from "react";

type PositionMode = "fixed" | "camera";

interface LightPositionData {
	camera_attached: boolean;
	x: number;
	y: number;
	z: number;
}

function isLightPosition(value: unknown): value is LightPositionData {
	return (
		typeof value === "object" &&
		value !== null &&
		"camera_attached" in value &&
		"x" in value &&
		"y" in value &&
		"z" in value
	);
}

const LightPositionRenderer = ({
	data,
	handleChange,
	path,
	label,
	required,
	errors,
	visible,
}: ControlProps) => {
	const posData = isLightPosition(data)
		? data
		: { camera_attached: true, x: 5, y: 2, z: 8 };

	const [mode, setMode] = useState<PositionMode>(
		posData.camera_attached ? "camera" : "fixed",
	);
	const [coords, setCoords] = useState<[string, string, string]>([
		String(posData.x),
		String(posData.y),
		String(posData.z),
	]);

	// Sync local state when external data changes
	useEffect(() => {
		if (!isLightPosition(data)) return;
		setMode(data.camera_attached ? "camera" : "fixed");
		setCoords([String(data.x), String(data.y), String(data.z)]);
	}, [data]);

	const emitChange = useCallback(
		(cameraAttached: boolean, x: number, y: number, z: number) => {
			handleChange(path, { camera_attached: cameraAttached, x, y, z });
		},
		[handleChange, path],
	);

	const handleModeChange = useCallback(
		(_: React.MouseEvent<HTMLElement>, newMode: PositionMode | null) => {
			if (newMode === null) return;
			setMode(newMode);

			const x = Number.parseFloat(coords[0]) || 0;
			const y = Number.parseFloat(coords[1]) || 0;
			const z = Number.parseFloat(coords[2]) || 0;
			emitChange(newMode === "camera", x, y, z);
		},
		[coords, emitChange],
	);

	const handleCoordChange = useCallback(
		(index: number, value: string) => {
			const newCoords = [...coords] as [string, string, string];
			newCoords[index] = value;
			setCoords(newCoords);

			const x = Number.parseFloat(newCoords[0]);
			const y = Number.parseFloat(newCoords[1]);
			const z = Number.parseFloat(newCoords[2]);

			if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
				emitChange(mode === "camera", x, y, z);
			}
		},
		[coords, mode, emitChange],
	);

	if (!visible) return null;

	return (
		<Box sx={{ marginBottom: 2 }}>
			<Typography
				variant="subtitle2"
				sx={{ marginBottom: 1, color: "text.secondary" }}
			>
				{label}
				{required && " *"}
			</Typography>

			<ToggleButtonGroup
				value={mode}
				exclusive
				onChange={handleModeChange}
				size="small"
				sx={{ marginBottom: 1.5 }}
			>
				<ToggleButton value="fixed">
					<LocationOnIcon sx={{ mr: 0.5 }} fontSize="small" />
					Fixed
				</ToggleButton>
				<ToggleButton value="camera">
					<CameraAltIcon sx={{ mr: 0.5 }} fontSize="small" />
					Camera
				</ToggleButton>
			</ToggleButtonGroup>

			<Paper elevation={1} sx={{ padding: 2 }}>
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
				<Typography
					variant="caption"
					color="text.secondary"
					sx={{ mt: 1, display: "block" }}
				>
					{mode === "fixed"
						? "Position in world coordinates"
						: "Offset from camera position"}
				</Typography>
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
 * Matches schemas that look like a LightPosition object:
 * has "camera_attached", "x", "y", "z" properties.
 */
export const lightPositionTester = rankWith(
	10,
	and(
		uiTypeIs("Control"),
		schemaMatches((schema) => {
			const props = (schema as any).properties;
			return (
				props &&
				"camera_attached" in props &&
				"x" in props &&
				"y" in props &&
				"z" in props
			);
		}),
	),
);

export default withJsonFormsControlProps(LightPositionRenderer);
