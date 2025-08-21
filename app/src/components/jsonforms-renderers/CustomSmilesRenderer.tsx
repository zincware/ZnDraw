import { useState } from "react";
import { withJsonFormsControlProps } from "@jsonforms/react";
import { FormLabel, Box, Tooltip, IconButton, TextField } from "@mui/material";
import { HelpOutline, Edit } from "@mui/icons-material";
import { SmilesModal } from "./SmilesModal";

const CustomSmilesRenderer = ({
	data,
	handleChange,
	path,
	label,
	schema,
	required,
	errors,
}: any) => {
	const [isModalOpen, setIsModalOpen] = useState(false);
	const value = data ?? schema.default ?? "";
	const hasError = errors && errors.length > 0;

	const handleSmilesChange = (newSmiles: string) => {
		handleChange(path, newSmiles);
		setIsModalOpen(false);
	};

	return (
		<Box sx={{ mb: 2 }}>
			<Box sx={{ display: "flex", alignItems: "center", mb: 1 }}>
				<FormLabel required={required}>{label}</FormLabel>
				{schema.description && (
					<Tooltip title={schema.description}>
						<IconButton size="small" sx={{ p: 0.5 }}>
							<HelpOutline fontSize="small" />
						</IconButton>
					</Tooltip>
				)}
			</Box>

			<Box sx={{ display: "flex", gap: 1, alignItems: "center" }}>
				<TextField
					value={value}
					onChange={(e) => handleChange(path, e.target.value)}
					placeholder="Enter SMILES string..."
					fullWidth
					variant="outlined"
					size="small"
					error={hasError}
					helperText={hasError ? errors[0].message : ""}
				/>
				<Tooltip title="Open Molecule Editor">
					<IconButton onClick={() => setIsModalOpen(true)} color="primary">
						<Edit />
					</IconButton>
				</Tooltip>
			</Box>

			{isModalOpen && (
				<SmilesModal
					smiles={value}
					setSmiles={handleSmilesChange}
					onClose={() => setIsModalOpen(false)}
				/>
			)}
		</Box>
	);
};

export default withJsonFormsControlProps(CustomSmilesRenderer);
