import FormControl from "@mui/material/FormControl";
import InputLabel from "@mui/material/InputLabel";
import MenuItem from "@mui/material/MenuItem";
import Select, { type SelectChangeEvent } from "@mui/material/Select";
import type { Visibility } from "../myapi/client";

interface Props {
	value: Visibility;
	onChange: (v: Visibility) => void;
	disabled?: boolean;
	size?: "small" | "medium";
	label?: string;
}

export default function VisibilitySelector({
	value,
	onChange,
	disabled,
	size = "small",
	label = "Visibility",
}: Props) {
	return (
		<FormControl size={size} disabled={disabled} fullWidth>
			<InputLabel>{label}</InputLabel>
			<Select
				value={value}
				label={label}
				onChange={(e: SelectChangeEvent) =>
					onChange(e.target.value as Visibility)
				}
			>
				<MenuItem value="private">Private — only me</MenuItem>
				<MenuItem value="group">Group — owning group members</MenuItem>
				<MenuItem value="public">Public — everyone</MenuItem>
			</Select>
		</FormControl>
	);
}
