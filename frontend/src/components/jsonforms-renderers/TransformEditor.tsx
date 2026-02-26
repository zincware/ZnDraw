import CloseIcon from "@mui/icons-material/Close";
import FilterAltIcon from "@mui/icons-material/FilterAlt";
import {
	Box,
	CircularProgress,
	FormControl,
	IconButton,
	InputLabel,
	MenuItem,
	Paper,
	Select,
	TextField,
	Tooltip,
	Typography,
} from "@mui/material";
import { useQuery } from "@tanstack/react-query";
import { getFrames } from "../../myapi/client";
import { useAppStore } from "../../store";
import type { Transform } from "../../utils/transformProcessor";

interface TransformEditorProps {
	value: Transform;
	label: string;
	required?: boolean;
	onChange: (newValue: Transform) => void;
	onClear: () => void;
}

interface ConstraintEntry {
	index: number;
	name: string;
	kwargs: Record<string, any>;
}

/**
 * Parse constraint data from frame into selectable entries.
 */
function parseConstraints(data: any): ConstraintEntry[] {
	if (!Array.isArray(data)) return [];
	return data
		.map((item: any, index: number) => {
			if (
				typeof item === "object" &&
				item !== null &&
				typeof item.name === "string" &&
				typeof item.kwargs === "object"
			) {
				return { index, name: item.name, kwargs: item.kwargs };
			}
			return null;
		})
		.filter(
			(entry: ConstraintEntry | null): entry is ConstraintEntry =>
				entry !== null,
		);
}

/**
 * Format a constraint entry for display in the selector.
 */
function formatConstraintLabel(entry: ConstraintEntry): string {
	const indices = entry.kwargs.indices;
	if (Array.isArray(indices)) {
		const preview =
			indices.length > 5
				? `[${indices.slice(0, 5).join(", ")}, ...]`
				: `[${indices.join(", ")}]`;
		return `#${entry.index}: ${entry.name} — indices: ${preview}`;
	}
	return `#${entry.index}: ${entry.name}`;
}

/**
 * TransformEditor component - inline editor for transform objects.
 *
 * When source is "constraints", provides a dropdown to select from
 * actual constraint entries in the current frame. Otherwise shows
 * text fields for source, path, and filter.
 */
export default function TransformEditor({
	value,
	label,
	required,
	onChange,
	onClear,
}: TransformEditorProps) {
	const transform = value || {
		type: "in_array" as const,
		source: "",
		path: "",
		filter: "",
	};

	const roomId = useAppStore((state) => state.roomId);
	const currentFrame = useAppStore((state) => state.currentFrame);

	// Fetch constraint data when source is "constraints"
	const isConstraintSource = transform.source === "constraints";
	const { data: constraintData, isLoading: isLoadingConstraints } = useQuery({
		queryKey: ["frame", roomId, currentFrame, "constraints"],
		queryFn: ({ signal }) =>
			getFrames(roomId!, currentFrame, ["constraints"], signal),
		enabled: isConstraintSource && !!roomId,
	});

	const constraints = isConstraintSource
		? parseConstraints(constraintData?.constraints)
		: [];

	// Extract current selection index from path (e.g., "0.kwargs.indices" -> 0)
	const selectedIndex = (() => {
		const match = transform.path.match(/^(\d+)\.kwargs\.indices$/);
		return match ? Number.parseInt(match[1], 10) : -1;
	})();

	const updateField = (field: string, newValue: string) => {
		onChange({
			...transform,
			[field]: newValue,
		} as Transform);
	};

	const handleConstraintSelect = (index: number) => {
		onChange({
			...transform,
			path: `${index}.kwargs.indices`,
		} as Transform);
	};

	return (
		<Box sx={{ marginBottom: 2 }}>
			<Box
				sx={{
					display: "flex",
					alignItems: "center",
					justifyContent: "space-between",
					marginBottom: 1,
				}}
			>
				<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
					<FilterAltIcon fontSize="small" color="primary" />
					<Typography variant="subtitle2" color="primary">
						{label}
						{required && " *"} (Transform)
					</Typography>
				</Box>
				<Tooltip title="Remove transform and use simple key">
					<IconButton size="small" onClick={onClear}>
						<CloseIcon fontSize="small" />
					</IconButton>
				</Tooltip>
			</Box>

			<Paper elevation={1} sx={{ padding: 2, backgroundColor: "#f9f9f9" }}>
				<Box sx={{ display: "flex", flexDirection: "column", gap: 1.5 }}>
					<TextField
						fullWidth
						size="small"
						label="Source"
						value={transform.source || ""}
						onChange={(e) => updateField("source", e.target.value)}
						placeholder="e.g., constraints"
						helperText="Frame data key containing indices"
					/>

					{isConstraintSource ? (
						<FormControl fullWidth size="small">
							<InputLabel>Constraint</InputLabel>
							<Select
								value={selectedIndex >= 0 ? selectedIndex : ""}
								label="Constraint"
								onChange={(e) =>
									handleConstraintSelect(e.target.value as number)
								}
								disabled={isLoadingConstraints}
								endAdornment={
									isLoadingConstraints ? (
										<CircularProgress size={20} sx={{ mr: 2 }} />
									) : undefined
								}
							>
								{constraints.length === 0 && !isLoadingConstraints ? (
									<MenuItem disabled>
										No constraints in current frame
									</MenuItem>
								) : (
									constraints.map((entry) => (
										<MenuItem key={entry.index} value={entry.index}>
											{formatConstraintLabel(entry)}
										</MenuItem>
									))
								)}
							</Select>
						</FormControl>
					) : (
						<TextField
							fullWidth
							size="small"
							label="Path"
							value={transform.path || ""}
							onChange={(e) => updateField("path", e.target.value)}
							placeholder="e.g., 0.kwargs.indices"
							helperText="Dot-separated path to extract indices"
						/>
					)}

					<TextField
						fullWidth
						size="small"
						label="Filter"
						value={transform.filter || ""}
						onChange={(e) => updateField("filter", e.target.value)}
						placeholder="e.g., arrays.positions"
						helperText="Frame data key to filter"
					/>
				</Box>
			</Paper>
		</Box>
	);
}
