import {
	Card,
	CardContent,
	Typography,
	Box,
	IconButton,
	Divider,
} from "@mui/material";
import CloseIcon from "@mui/icons-material/Close";
import { Rnd } from "react-rnd";
import type { Frame } from "./particles";
import { useSlowFrame } from "../contexts/SlowFrameContext";
import { PlotDisplay, useArrayPlots } from "./plotting";

export const ParticleInfoOverlay = ({
	show,
	info,
	position,
	particleIndex,
}: {
	show: boolean;
	info: { [key: string]: any };
	position: { x: number; y: number };
	particleIndex?: number;
}) => {
	const { atomsArrays } = useSlowFrame();
	const arrayPlots = useArrayPlots(atomsArrays, particleIndex ?? -1);

	if (!show) {
		return null;
	}

	const hasPlots = Object.keys(arrayPlots).length > 0;

	return (
		<Card
			sx={{
				position: "absolute",
				top: position.y - 25,
				left: position.x - 50,
				zIndex: 1000,
				maxWidth: hasPlots ? "24rem" : "18rem",
			}}
		>
			<CardContent sx={{ padding: "8px !important" }}>
				{Object.entries(info).map(([key, value]) => (
					<Typography variant="body2" key={key}>
						<strong>{key}:</strong> {String(value)}
					</Typography>
				))}

				{hasPlots && (
					<>
						<Divider sx={{ my: 1 }} />
						<Box sx={{ display: "flex", flexDirection: "column", gap: 1 }}>
							{Object.entries(arrayPlots).map(([key, plotInfo]) => (
								<Box
									key={key}
									sx={{ display: "flex", flexDirection: "column", gap: 0.5 }}
								>
									<Typography variant="caption" sx={{ fontWeight: "bold" }}>
										{key}:
									</Typography>
									<PlotDisplay
										plotData={plotInfo.plotData}
										plotType={plotInfo.plotType}
										plotLayout={plotInfo.plotLayout}
										width={180}
										height={120}
									/>
								</Box>
							))}
						</Box>
					</>
				)}
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
	const { atomsInfo } = useSlowFrame();

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
					<Typography variant="body2">
						Particles: {frame.positions.length}
					</Typography>
					<Typography variant="body2">
						Bonds: {frame.connectivity?.length || 0}
					</Typography>
					{atomsInfo?.energy && (
						<Typography variant="body2">
							Energy: {atomsInfo.energy.toFixed(3)} eV
						</Typography>
					)}
				</CardContent>
			</Card>
		</Rnd>
	);
};
