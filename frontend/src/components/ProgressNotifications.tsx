import CloseIcon from "@mui/icons-material/Close";
import {
	Box,
	CircularProgress,
	Fade,
	IconButton,
	LinearProgress,
	Paper,
	Typography,
} from "@mui/material";
import { useState } from "react";
import { LAYOUT_CONSTANTS } from "../constants/layout";
import type { Progress } from "../store";
import { useAppStore } from "../store";

/**
 * Format seconds as MM:SS (matching tqdm elapsed display).
 */
function formatElapsed(seconds: number): string {
	const m = Math.floor(seconds / 60);
	const s = Math.floor(seconds % 60);
	return `${String(m).padStart(2, "0")}:${String(s).padStart(2, "0")}`;
}

/**
 * Format rate with appropriate precision.
 */
function formatRate(rate: number): string {
	if (rate >= 100) return rate.toFixed(0);
	if (rate >= 10) return rate.toFixed(1);
	return rate.toFixed(2);
}

/**
 * Build a tqdm-like status string from tracker fields.
 *
 * Known total:   42/100 frames [00:05, 8.40 frames/s]
 * Unknown total: 42 it [00:05, 8.40 it/s]
 */
function formatStatus(tracker: Progress): string {
	const { n, total, elapsed, unit } = tracker;
	const rate = elapsed > 0 ? n / elapsed : 0;
	const ratePart = elapsed > 0 ? `, ${formatRate(rate)} ${unit}/s` : "";

	if (total !== null) {
		return `${n}/${total} ${unit} [${formatElapsed(elapsed)}${ratePart}]`;
	}
	return `${n} ${unit} [${formatElapsed(elapsed)}${ratePart}]`;
}

export default function ProgressNotifications() {
	const progressTrackers = useAppStore((state) => state.progressTrackers);
	const [hiddenProgressIds, setHiddenProgressIds] = useState<Set<string>>(
		new Set(),
	);

	const progressList = Object.values(progressTrackers).filter(
		(tracker) => !hiddenProgressIds.has(tracker.progress_id),
	);

	const handleHideProgress = (progressId: string) => {
		setHiddenProgressIds((prev) => new Set(prev).add(progressId));
	};

	if (progressList.length === 0) {
		return null;
	}

	return (
		<Box
			sx={{
				position: "fixed",
				bottom: LAYOUT_CONSTANTS.PROGRESSBAR_HEIGHT + 24,
				right: 24,
				zIndex: (theme) => theme.zIndex.drawer + 2,
				display: "flex",
				flexDirection: "column",
				gap: 1,
				maxWidth: 400,
			}}
		>
			{progressList.map((tracker) => {
				const hasDeterminate = tracker.total !== null && tracker.total > 0;
				const percent = hasDeterminate
					? (tracker.n / tracker.total!) * 100
					: undefined;

				return (
					<Fade key={tracker.progress_id} in={true}>
						<Paper
							elevation={8}
							sx={{
								p: 2,
								display: "flex",
								flexDirection: "column",
								gap: 1.5,
								minWidth: 320,
								position: "relative",
							}}
						>
							<IconButton
								size="small"
								onClick={() => handleHideProgress(tracker.progress_id)}
								sx={{
									position: "absolute",
									top: 8,
									right: 8,
									opacity: 0.6,
									"&:hover": { opacity: 1 },
								}}
							>
								<CloseIcon fontSize="small" />
							</IconButton>

							<Box sx={{ pr: 4 }}>
								<Typography variant="body1" fontWeight={600}>
									{tracker.description}
								</Typography>
							</Box>

							{hasDeterminate ? (
								<Box
									sx={{ display: "flex", flexDirection: "column", gap: 0.5 }}
								>
									<LinearProgress
										variant="determinate"
										value={percent}
										sx={{ height: 6, borderRadius: 3 }}
									/>
									<Typography
										variant="caption"
										color="text.secondary"
										fontFamily="monospace"
									>
										{formatStatus(tracker)}
									</Typography>
								</Box>
							) : (
								<Box sx={{ display: "flex", alignItems: "center", gap: 2 }}>
									<CircularProgress size={24} />
									<Typography
										variant="caption"
										color="text.secondary"
										fontFamily="monospace"
									>
										{formatStatus(tracker)}
									</Typography>
								</Box>
							)}
						</Paper>
					</Fade>
				);
			})}
		</Box>
	);
}
