import { type ControlProps, rankWith, schemaMatches } from "@jsonforms/core";
import { withJsonFormsControlProps } from "@jsonforms/react";
import LockIcon from "@mui/icons-material/Lock";
import LockOpenIcon from "@mui/icons-material/LockOpen";
import PersonIcon from "@mui/icons-material/Person";
import { Box, Button, FormLabel, Typography } from "@mui/material";
import { useAppStore } from "../../store";

const OwnershipToggle = ({ data, handleChange, path, label }: ControlProps) => {
	const userId = useAppStore((state) => state.user?.id ?? null);
	const isSuperuser = useAppStore((state) => state.user?.is_superuser ?? false);

	const isOwned = data != null && data !== "";
	const isOwnedByMe = isOwned && data === userId;
	const isOwnedByOther = isOwned && data !== userId;

	const handleClaim = () => {
		if (userId) {
			handleChange(path, userId);
		}
	};

	const handleRelease = () => {
		handleChange(path, null);
	};

	return (
		<Box sx={{ mb: 2 }}>
			<FormLabel>{label}</FormLabel>
			<Box sx={{ display: "flex", alignItems: "center", gap: 1, mt: 0.5 }}>
				{!isOwned && (
					<Button
						variant="outlined"
						size="small"
						startIcon={<PersonIcon />}
						onClick={handleClaim}
						disabled={!userId}
					>
						Claim ownership
					</Button>
				)}
				{isOwnedByMe && (
					<Button
						variant="outlined"
						size="small"
						color="warning"
						startIcon={<PersonIcon />}
						onClick={handleRelease}
					>
						Release ownership
					</Button>
				)}
				{isOwnedByOther && !isSuperuser && (
					<Typography
						variant="body2"
						color="text.secondary"
						sx={{ display: "flex", alignItems: "center", gap: 0.5 }}
					>
						<LockIcon fontSize="small" />
						Owned by another user
					</Typography>
				)}
				{isOwnedByOther && isSuperuser && (
					<>
						<Button
							variant="outlined"
							size="small"
							startIcon={<PersonIcon />}
							onClick={handleClaim}
						>
							Claim ownership
						</Button>
						<Button
							variant="outlined"
							size="small"
							color="warning"
							startIcon={<LockOpenIcon />}
							onClick={handleRelease}
						>
							Release ownership
						</Button>
					</>
				)}
			</Box>
		</Box>
	);
};

export const ownershipToggleTester = rankWith(
	10,
	schemaMatches(
		(schema) =>
			(schema as Record<string, unknown>)?.["x-custom-type"] ===
			"ownership-toggle",
	),
);

export default withJsonFormsControlProps(OwnershipToggle);
