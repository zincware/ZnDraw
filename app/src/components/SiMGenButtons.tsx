import { useEffect, useState } from "react";
import Box from "@mui/material/Box";
import Tooltip from "@mui/material/Tooltip";
import Button from "@mui/material/Button";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";
import LinkIcon from "@mui/icons-material/Link";
import ScienceIcon from "@mui/icons-material/Science";
import ArticleIcon from "@mui/icons-material/Article";
import { submitExtension } from "../myapi/client";
import { useAppStore } from "../store";

interface SiMGenButtonsProps {
	onTutorialClick: () => void;
}

export default function SiMGenButtons({ onTutorialClick }: SiMGenButtonsProps) {
	// Use individual selectors to prevent unnecessary re-renders
	const roomId = useAppStore((state) => state.roomId);
	const globalSettings = useAppStore((state) => state.globalSettings);
	const showSnackbar = useAppStore((state) => state.showSnackbar);
	const [showAnimation, setShowAnimation] = useState(true);

	// Disable animation after it plays once
	useEffect(() => {
		if (showAnimation) {
			const timer = setTimeout(() => {
				setShowAnimation(false);
			}, 3000); // Animation duration (3 pulses of 1s each)
			return () => clearTimeout(timer);
		}
	}, [showAnimation]);

	// Handler for submitting SiMGen extensions
	const handleExtensionSubmit = async (extensionName: string) => {
		if (!roomId) {
			showSnackbar("Please select a room first", "warning");
			return;
		}

		try {
			// SiMGen extensions are registered publicly
			await submitExtension({
				roomId,
				category: "modifiers",
				extension: extensionName,
				data: {},
				isPublic: true,
			});
			showSnackbar(`${extensionName} submitted successfully`, "success");
		} catch (error) {
			console.error(`Failed to submit ${extensionName}:`, error);
			showSnackbar(`Failed to submit ${extensionName}`, "error");
		}
	};

	// Don't render anything if simgen is not enabled
	if (!globalSettings?.simgen?.enabled) {
		return null;
	}

	return (
		<Box
			sx={{
				display: "flex",
				alignItems: "center",
				gap: 1.5,
				m: 0,
			}}
		>
			<Tooltip title="Connect selected atoms">
				<Button
					color="inherit"
					startIcon={<LinkIcon />}
					onClick={() => handleExtensionSubmit("Connect")}
					size="medium"
					sx={{
						textTransform: "none",
						fontWeight: 500,
						"&:hover": {
							bgcolor: (theme) =>
								theme.palette.mode === "light"
									? "rgba(25, 118, 210, 0.08)"
									: "rgba(144, 202, 249, 0.12)",
						},
					}}
				>
					Connect
				</Button>
			</Tooltip>
			<Tooltip title="Run SiMGen Demo">
				<Button
					color="inherit"
					startIcon={<ScienceIcon />}
					onClick={() => handleExtensionSubmit("SiMGenDemo")}
					size="medium"
					sx={{
						textTransform: "none",
						fontWeight: 500,
						"&:hover": {
							bgcolor: (theme) =>
								theme.palette.mode === "light"
									? "rgba(25, 118, 210, 0.08)"
									: "rgba(144, 202, 249, 0.12)",
						},
					}}
				>
					Generate
				</Button>
			</Tooltip>
			<Tooltip title="Replace scene with empty canvas">
				<Button
					color="inherit"
					startIcon={<ArticleIcon />}
					onClick={() => handleExtensionSubmit("NewCanvas")}
					size="medium"
					sx={{
						textTransform: "none",
						fontWeight: 500,
						"&:hover": {
							bgcolor: (theme) =>
								theme.palette.mode === "light"
									? "rgba(25, 118, 210, 0.08)"
									: "rgba(144, 202, 249, 0.12)",
						},
					}}
				>
					New Canvas
				</Button>
			</Tooltip>
			<Box
				sx={{
					width: "1px",
					height: "24px",
					bgcolor: "divider",
					mx: 0.5,
				}}
			/>
			<Tooltip title="SiMGen Tutorial">
				<Button
					color="inherit"
					onClick={onTutorialClick}
					size="medium"
					sx={{
						minWidth: "auto",
						px: 1,
						textTransform: "none",
						fontWeight: 500,
						"&:hover": {
							bgcolor: (theme) =>
								theme.palette.mode === "light"
									? "rgba(25, 118, 210, 0.08)"
									: "rgba(144, 202, 249, 0.12)",
						},
					}}
				>
					<HelpOutlineIcon fontSize="small" />
				</Button>
			</Tooltip>
		</Box>
	);
}
