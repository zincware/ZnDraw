import React, { useState, useEffect, useMemo } from "react";
import { useAppContext } from "../../contexts/AppContext";
import HeadBar from "../headbar";
import Sidebar from "../sidebar";
import FrameProgressBar from "../progressbar";
import { Plotting } from "../plotting";
import { ParticleInfoOverlay, SceneInfoOverlay } from "../overlays";
import * as THREE from "three";

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

		// Selection and interaction
		selectedIds,
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
			{/* Main UI Components */}
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

			<Sidebar token={token} />

			<FrameProgressBar
				step={step}
				setStep={setStep}
				length={length}
				selectedFrames={selectedFrames}
				setSelectedFrames={setSelectedFrames}
				bookmarks={bookmarks}
				setBookmarks={setBookmarks}
				connected={connected}
			/>

			{/* Plotting Component */}
			{updatedPlotsList.length > 0 && (
				<Plotting token={token} updatedPlotsList={updatedPlotsList} />
			)}

			{/* Overlays */}
			{showParticleInfo &&
				hoveredId !== -1 &&
				hoveredId < currentFrame.positions.length && (
					<ParticleInfoOverlay
						particle={{
							position: currentFrame.positions[hoveredId],
							number: currentFrame.numbers[hoveredId],
							color: currentFrame.arrays.colors[hoveredId],
							radius: currentFrame.arrays.radii[hoveredId],
							...(isDrawing && { Line: `${lineLength.toFixed(2)} Å` }),
						}}
						position={cursorPosition}
					/>
				)}

			<SceneInfoOverlay
				frame={currentFrame}
				setShowParticleInfo={setShowParticleInfo}
			/>
		</>
	);
};
