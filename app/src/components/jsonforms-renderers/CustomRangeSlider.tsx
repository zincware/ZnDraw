// src/renderers/CustomRangeSlider.js
import { withJsonFormsControlProps } from "@jsonforms/react";
import {
	FormLabel,
	Slider,
	Box,
	Typography,
	Tooltip,
	IconButton,
} from "@mui/material";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";

const CustomRangeSlider = ({
	data,
	handleChange,
	path,
	label,
	schema,
	required,
}) => {
	const step = schema.step || 0.1;
	const value = data ?? schema.default ?? schema.minimum;

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
			<Box sx={{ display: "flex", alignItems: "center", gap: 2, paddingX: 1 }}>
				<Slider
					value={value}
					onChange={(_event, newValue) => handleChange(path, newValue)}
					min={schema.minimum}
					max={schema.maximum}
					step={step}
					aria-labelledby="input-slider"
				/>
				<Typography sx={{ minWidth: 35 }}>{value}</Typography>
			</Box>
		</Box>
	);
};

export default withJsonFormsControlProps(CustomRangeSlider);
