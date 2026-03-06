import {
	Box,
	Chip,
	FormControl,
	InputLabel,
	MenuItem,
	Select,
	Typography,
} from "@mui/material";
import React from "react";
import type { ProviderInfo } from "../../myapi/client";

interface FilesystemSelectorProps {
	providers: ProviderInfo[];
	selected: string;
	onSelect: (fullName: string) => void;
}

/**
 * Dropdown selector for choosing a filesystem provider.
 * Uses provider full_name as the unique key.
 */
export function FilesystemSelector({
	providers,
	selected,
	onSelect,
}: FilesystemSelectorProps) {
	if (providers.length === 0) {
		return null;
	}

	// If only one provider, show info instead of dropdown
	if (providers.length === 1) {
		const p = providers[0];
		return (
			<Box sx={{ mb: 2 }}>
				<Typography variant="body2" color="text.secondary">
					<strong>Filesystem:</strong> {p.name}
				</Typography>
			</Box>
		);
	}

	return (
		<Box sx={{ mb: 3 }}>
			<FormControl fullWidth>
				<InputLabel>Select Filesystem</InputLabel>
				<Select
					value={selected}
					label="Select Filesystem"
					onChange={(e) => onSelect(e.target.value)}
				>
					{providers.map((p) => (
						<MenuItem key={p.id} value={p.full_name}>
							<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
								<Typography>{p.name}</Typography>
								<Chip
									label={p.category}
									size="small"
									color="primary"
									variant="outlined"
								/>
							</Box>
						</MenuItem>
					))}
				</Select>
			</FormControl>
		</Box>
	);
}
