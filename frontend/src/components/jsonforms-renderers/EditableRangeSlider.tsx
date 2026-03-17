import { type ControlProps, rankWith, schemaMatches } from "@jsonforms/core";
import {
	useJsonForms,
	withJsonFormsControlProps,
} from "@jsonforms/react";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";
import {
	Box,
	FormLabel,
	IconButton,
	Slider,
	TextField,
	Tooltip,
	Typography,
} from "@mui/material";
import { useCallback, useState } from "react";

/**
 * Editable-range slider: a slider whose min/max labels are click-to-edit.
 *
 * Schema annotations consumed:
 *   x-custom-type: "editable-range"
 *   x-min-field:   sibling field name for the minimum bound
 *   x-max-field:   sibling field name for the maximum bound
 *   step:          slider step size
 */
const EditableRangeSlider = ({
	data,
	handleChange,
	path,
	label,
	schema,
	required,
}: ControlProps) => {
	const ctx = useJsonForms();
	const rootData = ctx.core?.data ?? {};

	const ext = schema as Record<string, unknown>;
	const minField = ext["x-min-field"] as string;
	const maxField = ext["x-max-field"] as string;
	const step = (ext.step as number) || 0.001;

	const min = Number(rootData[minField] ?? -0.25);
	const max = Number(rootData[maxField] ?? 0.25);
	const value = Number(data ?? schema.default ?? min);

	const [editingMin, setEditingMin] = useState(false);
	const [editingMax, setEditingMax] = useState(false);
	const [editingValue, setEditingValue] = useState(false);
	const [minDraft, setMinDraft] = useState("");
	const [maxDraft, setMaxDraft] = useState("");
	const [valueDraft, setValueDraft] = useState("");

	const commitMin = useCallback(() => {
		const parsed = Number.parseFloat(minDraft);
		if (!Number.isNaN(parsed)) {
			handleChange(minField, parsed);
			if (value < parsed) {
				handleChange(path, parsed);
			}
		}
		setEditingMin(false);
	}, [minDraft, handleChange, minField, path, value]);

	const commitMax = useCallback(() => {
		const parsed = Number.parseFloat(maxDraft);
		if (!Number.isNaN(parsed)) {
			handleChange(maxField, parsed);
			if (value > parsed) {
				handleChange(path, parsed);
			}
		}
		setEditingMax(false);
	}, [maxDraft, handleChange, maxField, path, value]);

	const commitValue = useCallback(() => {
		const parsed = Number.parseFloat(valueDraft);
		if (!Number.isNaN(parsed)) {
			const clamped = Math.min(Math.max(parsed, min), max);
			handleChange(path, clamped);
		}
		setEditingValue(false);
	}, [valueDraft, handleChange, path, min, max]);

	return (
		<Box sx={{ marginBottom: 2 }}>
			<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
				<FormLabel required={required}>{label}</FormLabel>
				{schema.description && (
					<Tooltip title={schema.description} placement="top">
						<IconButton size="small" sx={{ padding: 0.5 }}>
							<HelpOutlineIcon fontSize="small" />
						</IconButton>
					</Tooltip>
				)}
			</Box>
			<Box sx={{ display: "flex", alignItems: "center", gap: 1, paddingX: 1 }}>
				{/* Min label — click to edit */}
				{editingMin ? (
					<TextField
						size="small"
						variant="standard"
						value={minDraft}
						onChange={(e) => setMinDraft(e.target.value)}
						onBlur={commitMin}
						onKeyDown={(e) => e.key === "Enter" && commitMin()}
						autoFocus
						sx={{ width: 64 }}
						slotProps={{ htmlInput: { style: { textAlign: "center" } } }}
					/>
				) : (
					<Typography
						sx={{
							minWidth: 40,
							cursor: "pointer",
							userSelect: "none",
							"&:hover": { color: "primary.main" },
						}}
						onClick={() => {
							setMinDraft(String(min));
							setEditingMin(true);
						}}
					>
						{min}
					</Typography>
				)}

				<Slider
					value={value}
					onChange={(_event, newValue) => handleChange(path, newValue)}
					min={min}
					max={max}
					step={step}
					sx={{ flexGrow: 1 }}
					aria-labelledby="editable-range-slider"
				/>

				{/* Max label — click to edit */}
				{editingMax ? (
					<TextField
						size="small"
						variant="standard"
						value={maxDraft}
						onChange={(e) => setMaxDraft(e.target.value)}
						onBlur={commitMax}
						onKeyDown={(e) => e.key === "Enter" && commitMax()}
						autoFocus
						sx={{ width: 64 }}
						slotProps={{ htmlInput: { style: { textAlign: "center" } } }}
					/>
				) : (
					<Typography
						sx={{
							minWidth: 40,
							cursor: "pointer",
							userSelect: "none",
							"&:hover": { color: "primary.main" },
						}}
						onClick={() => {
							setMaxDraft(String(max));
							setEditingMax(true);
						}}
					>
						{max}
					</Typography>
				)}

				{/* Value label — click to edit */}
				{editingValue ? (
					<TextField
						size="small"
						variant="standard"
						value={valueDraft}
						onChange={(e) => setValueDraft(e.target.value)}
						onBlur={commitValue}
						onKeyDown={(e) => e.key === "Enter" && commitValue()}
						autoFocus
						sx={{ width: 64 }}
						slotProps={{ htmlInput: { style: { textAlign: "center" } } }}
					/>
				) : (
					<Typography
						sx={{
							minWidth: 35,
							cursor: "pointer",
							userSelect: "none",
							"&:hover": { color: "primary.main" },
						}}
						onClick={() => {
							setValueDraft(String(value));
							setEditingValue(true);
						}}
					>
						{value}
					</Typography>
				)}
			</Box>
		</Box>
	);
};

export const editableRangeSliderTester = rankWith(
	10, // Higher than customRangeSlider (5) to take priority
	schemaMatches((schema) => {
		const ext = schema as Record<string, unknown>;
		return ext?.["x-custom-type"] === "editable-range";
	}),
);

export default withJsonFormsControlProps(EditableRangeSlider);
