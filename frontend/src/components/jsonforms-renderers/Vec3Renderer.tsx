/**
 * Vec3Renderer - Custom JSONForms renderer for vec3 fields.
 *
 * Matches schemas with `x-custom-type: "vec3"`. Handles both:
 * - Nested `[[x, y, z]]` (list[Vec3], e.g. CircleCurve)
 * - Flat `[x, y, z]` (tuple[float,float,float], e.g. Camera.up, Floor.position)
 */

import {
	type ControlProps,
	and,
	rankWith,
	schemaMatches,
	uiTypeIs,
} from "@jsonforms/core";
import { withJsonFormsControlProps } from "@jsonforms/react";
import { Box, TextField, Typography } from "@mui/material";
import { useCallback, useEffect, useState } from "react";

/** Whether data is nested `[[x, y, z]]` vs flat `[x, y, z]`. */
function isNested(data: unknown): boolean {
	return Array.isArray(data) && data.length === 1 && Array.isArray(data[0]);
}

function extractVec3(data: unknown): [number, number, number] {
	if (Array.isArray(data)) {
		if (isNested(data) && data[0].length === 3) {
			return data[0] as [number, number, number];
		}
		if (data.length === 3 && !Array.isArray(data[0])) {
			return data as [number, number, number];
		}
	}
	return [0, 0, 0];
}

const Vec3Renderer = ({
	data,
	handleChange,
	path,
	label,
	visible,
	errors,
}: ControlProps) => {
	const vec = extractVec3(data);

	const [coords, setCoords] = useState<[string, string, string]>([
		String(vec[0]),
		String(vec[1]),
		String(vec[2]),
	]);

	// Sync local state when external data changes
	useEffect(() => {
		const v = extractVec3(data);
		setCoords([String(v[0]), String(v[1]), String(v[2])]);
	}, [data]);

	const handleCoordChange = useCallback(
		(index: number, value: string) => {
			const newCoords = [...coords] as [string, string, string];
			newCoords[index] = value;
			setCoords(newCoords);

			const x = Number.parseFloat(newCoords[0]);
			const y = Number.parseFloat(newCoords[1]);
			const z = Number.parseFloat(newCoords[2]);

			if (!Number.isNaN(x) && !Number.isNaN(y) && !Number.isNaN(z)) {
				handleChange(path, isNested(data) ? [[x, y, z]] : [x, y, z]);
			}
		},
		[coords, data, handleChange, path],
	);

	if (!visible) return null;

	const labels = ["X", "Y", "Z"];

	return (
		<Box sx={{ marginBottom: 2 }}>
			<Typography
				variant="subtitle2"
				sx={{ marginBottom: 1, color: "text.secondary" }}
			>
				{label}
			</Typography>
			<Box sx={{ display: "flex", gap: 1 }}>
				{labels.map((axis, i) => (
					<TextField
						key={axis}
						label={axis}
						type="number"
						size="small"
						value={coords[i]}
						onChange={(e) => handleCoordChange(i, e.target.value)}
						error={!!errors}
						inputProps={{ step: "any" }}
						sx={{ flex: 1 }}
					/>
				))}
			</Box>
			{errors && (
				<Typography variant="caption" color="error" sx={{ mt: 0.5 }}>
					{errors}
				</Typography>
			)}
		</Box>
	);
};

export const vec3Tester = rankWith(
	10,
	and(
		uiTypeIs("Control"),
		schemaMatches(
			(schema) => (schema as Record<string, unknown>)["x-custom-type"] === "vec3",
		),
	),
);

export default withJsonFormsControlProps(Vec3Renderer);
