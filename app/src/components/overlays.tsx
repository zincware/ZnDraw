import { Card, CardContent, Typography, Box, IconButton } from "@mui/material";
import CloseIcon from "@mui/icons-material/Close";
import { Rnd } from "react-rnd";
import type { Frame } from "./particles";

export const ParticleInfoOverlay = ({
	show,
	info,
	position,
}: {
	show: boolean;
	info: { [key: string]: any };
	position: { x: number; y: number };
}) => {
	if (!show) {
		return null;
	}

	return (
		<Card
			sx={{
				position: "absolute",
				top: position.y - 25,
				left: position.x - 50,
				zIndex: 1000,
				maxWidth: "18rem",
			}}
		>
			<CardContent sx={{ padding: "8px !important" }}>
				{Object.entries(info).map(([key, value]) => (
					<Typography variant="body2" key={key}>
						<strong>{key}:</strong> {String(value)}
					</Typography>
				))}
			</CardContent>
		</Card>
	);
};

export const SceneInfoOverlay = ({
	frame,
	setShowParticleInfo,
}: {
	frame: Frame;
	setShowParticleInfo: (show: boolean) => void;
}) => {
	return (
		<Rnd
			default={{
				x: window.innerWidth / 2 - 150,
				y: 50,
				width: 280,
				height: 150,
			}}
			minHeight={150}
			minWidth={200}
			style={{ zIndex: 1000 }}
			dragHandleClassName="drag-handle" // Class for the drag handle
		>
			<Card
				sx={{
					width: "100%",
					height: "100%",
					display: "flex",
					flexDirection: "column",
				}}
			>
				{/* Header Box serves as the drag handle and contains title/button */}
				<Box
					className="drag-handle"
					sx={{
						display: "flex",
						justifyContent: "space-between",
						alignItems: "center",
						padding: "4px 8px",
						borderBottom: 1,
						borderColor: "divider",
						cursor: "move",
					}}
				>
					<Typography variant="h6" sx={{ fontSize: "1rem" }}>
						Info
					</Typography>
					<IconButton
						size="small"
						onClick={() => setShowParticleInfo(false)}
						aria-label="close"
					>
						<CloseIcon fontSize="small" />
					</IconButton>
				</Box>
				<CardContent
					sx={{
						textAlign: "start",
						overflowY: "auto",
						padding: "8px",
					}}
				>
					{frame.calc?.energy && (
						<Typography variant="body2">
							Energy: {frame.calc.energy.toFixed(3)} eV
						</Typography>
					)}
					<Typography variant="body2">
						Particles: {frame.positions.length}
					</Typography>
					<Typography variant="body2">
						Bonds: {frame.connectivity?.length || 0}
					</Typography>
				</CardContent>
			</Card>
		</Rnd>
	);
};
