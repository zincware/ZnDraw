import { Box } from "@mui/material";
import type React from "react";
import { useEffect, useMemo, useState } from "react";
import * as THREE from "three";
import { useAppContext } from "../../contexts/AppContext";
import HeadBar from "../headbar";
import { ParticleInfoOverlay, SceneInfoOverlay } from "../overlays";
import { Plotting } from "../plotting";
import FrameProgressBar from "../progressbar";
import Sidebar from "../sidebar";

export const UIContainer: React.FC = () => {
	const {
		// Connection and room state
		connected,
		token,
		roomName,
		roomLock,

		// Frame and player state
		currentFrame,
		step,
		setStep,
		length,
		playing,
		setPlaying,
		frameRate,
		setFrameRate,
		isFrameRendering,

		// Selection and interaction
		selectedIds,
		setSelectedIds,
		selectedFrames,
		setSelectedFrames,
		hoveredId,
		isDrawing,
		setIsDrawing,
		points,
		setPoints,
		selectedPoint,
		editMode,

		// UI state
		colorMode,
		handleColorMode,
		showParticleInfo,
		setShowParticleInfo,
		tutorialURL,
		showSiMGen,

		// Additional UI state
		modifierQueue,
		isAuthenticated,
		addPlotsWindow,
		setAddPlotsWindow,

		// Data collections
		bookmarks,
		setBookmarks,
		geometries,
		setGeometries,

		// Configuration
		roomConfig,

		// Plotting
		updatedPlotsList,

		// Messages
		messages,
		setMessages,
	} = useAppContext();

	// Local UI state
	const [cursorPosition, setCursorPosition] = useState({ x: 0, y: 0 });

	// Update cursor position for overlays
	useEffect(() => {
		const updateCursorPosition = (event: MouseEvent) => {
			setCursorPosition({ x: event.clientX, y: event.clientY });
		};

		window.addEventListener("mousemove", updateCursorPosition);

		return () => {
			window.removeEventListener("mousemove", updateCursorPosition);
		};
	}, []);

	// Calculate line length for drawing overlay
	const lineLength = useMemo(() => {
		if (points.length >= 2) {
			const start = points[points.length - 2];
			const end = points[points.length - 1];
			return start.distanceTo(end);
		}
		return 0;
	}, [points]);

	return (
		<>
			<Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
				<HeadBar
					room={roomName}
					colorMode={colorMode}
					handleColorMode={handleColorMode}
					setIsDrawing={setIsDrawing}
					setGeometries={setGeometries}
					setPoints={setPoints}
					isDrawing={isDrawing}
					tutorialURL={tutorialURL}
					showSiMGen={showSiMGen}
					modifierQueue={modifierQueue}
					isAuthenticated={isAuthenticated}
					roomLock={roomLock}
					setAddPlotsWindow={setAddPlotsWindow}
					messages={messages}
					setMessages={setMessages}
					token={token}
					step={step}
					selection={selectedIds}
				/>
				<Box sx={{ display: "flex", flexGrow: 1 }}>
					<Sidebar token={token} />
					<Box sx={{ flexGrow: 1 }} />
				</Box>
				<Box sx={{
                position: 'fixed', // Pins the component to the viewport
                bottom: 0,          // Aligns to the bottom
                left: 0,            // Aligns to the left
                right: 0,           // Extends to the right, ensuring full width
                zIndex: 1000,       // Ensure it's on top of other content
            }}>	
				<FrameProgressBar
					step={step}
					setStep={setStep}
					length={length}
					selectedFrames={selectedFrames}
					setSelectedFrames={setSelectedFrames}
					bookmarks={bookmarks}
					setBookmarks={setBookmarks}
					connected={connected}
					frameRate={frameRate}
					setFrameRate={setFrameRate}
					isFrameRendering={isFrameRendering}
				/>
				</Box>
			</Box>

			{/* Plotting Component */}
			<Plotting
				token={token}
				updatedPlotsList={updatedPlotsList}
				addPlotsWindow={addPlotsWindow}
				setStep={setStep}
				step={step}
				setSelectedFrames={setSelectedFrames}
				setSelectedIds={setSelectedIds}
			/>

			{/* Overlays */}
			{(() => {
				return (
					showParticleInfo &&
					hoveredId !== -1 &&
					hoveredId < currentFrame.positions.length && (
						<ParticleInfoOverlay
							show={showParticleInfo}
							info={{
								position: currentFrame.positions[hoveredId],
								number: currentFrame.numbers[hoveredId],
								color: currentFrame.arrays.colors[hoveredId],
								radius: currentFrame.arrays.radii[hoveredId],
								...(isDrawing && { Line: `${lineLength.toFixed(2)} Å` }),
							}}
							position={cursorPosition}
						/>
					)
				);
			})()}
			{showParticleInfo && (
				<SceneInfoOverlay
					frame={currentFrame}
					setShowParticleInfo={setShowParticleInfo}
				/>
			)}
		</>
	);
};
