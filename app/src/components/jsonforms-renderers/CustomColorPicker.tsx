// src/renderers/CustomColorPicker.js
import React from "react";
import { withJsonFormsControlProps } from "@jsonforms/react";
import { FormLabel, Box, Tooltip, IconButton } from "@mui/material";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";

const CustomColorPicker = ({
	data,
	handleChange,
	path,
	label,
	schema,
	required,
}) => {
	const value = data ?? schema.default ?? "#000000";

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
			<Box sx={{ display: "flex", alignItems: "center", gap: 2, marginTop: 1 }}>
				<input
					type="color"
					value={value}
					onChange={(e) => handleChange(path, e.target.value)}
					style={{
						width: "40px",
						height: "40px",
						border: "none",
						background: "none",
						cursor: "pointer",
					}}
				/>
				<code>{value}</code>
			</Box>
		</Box>
	);
};

export default withJsonFormsControlProps(CustomColorPicker);
