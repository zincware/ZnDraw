import { TextField, Tooltip } from "@mui/material";
import React, { useEffect, useState } from "react";
import { IoMdCodeDownload } from "react-icons/io";

interface FrameRateControlProps {
	frameRate: number;
	setFrameRate: (val: number) => void;
	isFrameRendering: boolean;
	connected: boolean;
}

export const FrameRateControl: React.FC<FrameRateControlProps> = ({
	frameRate,
	setFrameRate,
	isFrameRendering,
	connected,
}) => {
	const [showLoadingIcon, setShowLoadingIcon] = useState(false);

	// Delay before showing ⧗ to avoid flickering
	useEffect(() => {
		let timeout: NodeJS.Timeout;
		if (isFrameRendering) {
			timeout = setTimeout(() => setShowLoadingIcon(true), 50);
		} else {
			setShowLoadingIcon(false);
			clearTimeout(timeout);
		}
		return () => clearTimeout(timeout);
	}, [isFrameRendering]);

	const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
		const val = Number.parseInt(e.target.value, 10);
		if (!isNaN(val)) {
			setFrameRate(Math.max(1, val)); // Ensure frame rate is at least 1
		}
	};

	return (
		<div
			className="d-flex justify-content-center align-items-center user-select-none"
			style={{ display: "flex", height: "100%" }}
		>
			<TextField
				type="number"
				value={frameRate}
				onChange={handleChange}
				InputProps={{
					inputProps: {
						min: 1,
						style: {
							textAlign: "center",
							fontSize: "10px",
						},
					},
				}}
				size="small"
				style={{
					width: "40px",
					height: "18px",
				}}
				title="Frame rate (1 = every frame, 2 = every 2nd frame, etc.)"
			/>

			{/* Frame loading icon */}
			<Tooltip title="Loading frame data...">
				<div
					style={{
						height: "100%",
						animation: showLoadingIcon ? "pulse 1s infinite" : "none",
						visibility: showLoadingIcon ? "visible" : "hidden",
						alignItems: "center",
						justifyContent: "center",
						display: "flex",
					}}
					title="Loading frame data..."
				>
					<IoMdCodeDownload />
				</div>
			</Tooltip>

			{/* Server connection spinner */}
			{!connected && (
				<Tooltip title="Not connected to server">
					<div
						style={{
							height: "25px",
							width: "25px",
						}}
					>
						<div
							className="spinner-border spinner-border-sm text-primary"
							role="status"
							style={{
								width: "15px",
								height: "15px",
							}}
						></div>
						<span className="visually-hidden">Loading...</span>
					</div>
				</Tooltip>
			)}
		</div>
	);
};