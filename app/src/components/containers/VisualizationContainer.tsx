import React, { useMemo, useEffect } from "react";
import { Canvas } from "@react-three/fiber";
import { Pathtracer } from "@react-three/gpu-pathtracer";
import { Environment } from "@react-three/drei";
import { useAppContext } from "../../contexts/AppContext";
import {
	BondInstances,
	ParticleInstances,
	PerParticleVectors,
	SimulationCell,
	Player,
} from "../particles";
import { Geometries } from "../geometries";
import CameraAndControls from "../cameraAndControls";
import { Floor } from "../floor";
import { Line3D, VirtualCanvas } from "../lines";
import ControlsBuilder from "../transforms";
import VectorField from "../vectorfield";

interface VisualizationContainerProps {
	onPointerMissed: () => void;
	onDragOver: (event: React.DragEvent) => void;
	onDrop: (event: React.DragEvent) => void;
}

export const VisualizationContainer: React.FC<VisualizationContainerProps> = ({
	onPointerMissed,
	onDragOver,
	onDrop,
}) => {
	const {
		// Frame and player state
		currentFrame,
		setCurrentFrame,
		step,
		setStep,
		length,
		playing,
		setPlaying,
		isFrameRendering,
		setIsFrameRendering,
		frameRate,

		// Selection and interaction
		selectedIds,
		setSelectedIds,
		selectedFrames,
		hoveredId,
		setHoveredId,
		isDrawing,
		setIsDrawing,
		points,
		setPoints,
		selectedPoint,
		setSelectedPoint,
		editMode,
		setEditMode,

		// Camera and visual settings
		cameraAndControls,
		setCameraAndControls,
		colorMode,
		roomConfig,

		// Room state
		token,
		roomLock,

		// Geometry data
		geometries,

		// Vector data
		perParticleVectors,
		vectorFieldData,
		vectorColormap,
	} = useAppContext();

	// Track roomConfig changes for real-time updates
	useEffect(() => {
		console.log("UI updated with new config");
	}, [roomConfig]);

	// Extract constrained atom indices from frame constraints
	const constrainedAtomIds = useMemo(() => {
		const constrainedIds = new Set<number>();

		if (currentFrame?.constraints) {
			for (const constraint of currentFrame.constraints) {
				if (constraint.type === "FixAtoms" && constraint.indices) {
					for (const index of constraint.indices) {
						constrainedIds.add(index);
					}
				}
			}
		}

		return constrainedIds;
	}, [currentFrame?.constraints]);

	// Note: Frame rendering state is now tracked by websocket data loading in setupFrames

	return (
		<div
			className="canvas-container edit-mode-border"
			onDragOver={onDragOver}
			onDrop={onDrop}
			tabIndex={0}
		>
			<Canvas onPointerMissed={onPointerMissed} shadows>
				<CameraAndControls
					cameraConfig={roomConfig.Camera}
					cameraAndControls={cameraAndControls}
					setCameraAndControls={setCameraAndControls}
					currentFrame={currentFrame}
					selectedIds={selectedIds}
					colorMode={colorMode}
				/>

				<Pathtracer enabled={roomConfig.PathTracer.enabled}>
					{roomConfig.PathTracer.enabled &&
						roomConfig.PathTracer.environment !== "none" && (
							<Environment preset={roomConfig.PathTracer.environment as any} />
						)}

					{/* Lighting */}
					{roomConfig.Visualization.floor ? (
						<>
							<Floor colorMode={colorMode} roomConfig={roomConfig} />
							<directionalLight
								position={[0, 100, 0]}
								intensity={1.0}
								castShadow
								shadow-mapSize-width={
									(roomConfig.Camera.camera_far || 300) * 10
								} // Adjust the width of the shadow map
								shadow-mapSize-height={
									(roomConfig.Camera.camera_far || 300) * 10
								} // Adjust the height of the shadow map
								shadow-camera-near={10} // Adjust the near clipping plane of the shadow camera
								shadow-camera-far={800} // Adjust the far clipping plane of the shadow camera
								shadow-camera-left={-1 * (roomConfig.Camera.camera_far || 300)} // Set the left boundary for the shadow camera frustum
								shadow-camera-right={roomConfig.Camera.camera_far || 300} // Set the right boundary for the shadow camera frustum
								shadow-camera-top={roomConfig.Camera.camera_far || 300} // Set the top boundary for the shadow camera frustum
								shadow-camera-bottom={
									-1 * (roomConfig.Camera.camera_far || 300)
								} // Set the bottom boundary for the shadow camera frustum
							/>
						</>
					) : (
						<directionalLight position={[0, 100, 0]} intensity={1.0} />
					)}

					{/* Vector field */}
					{roomConfig.VectorDisplay.vectorfield &&
						vectorFieldData.length > 0 && (
							<VectorField
								vectors={vectorFieldData}
								arrowsConfig={{
									normalize: roomConfig.VectorDisplay.normalize || true,
									colorrange: (roomConfig.VectorDisplay.colorrange || [
										0, 1.0,
									]) as [number, number],
									scale_vector_thickness:
										roomConfig.VectorDisplay.scale_vector_thickness || false,
									opacity: roomConfig.VectorDisplay.opacity || 1.0,
									colormap:
										vectorColormap.length > 0
											? vectorColormap
											: (roomConfig.VectorDisplay.default_colormap as [
													number,
													number,
													number,
												][]) || [
													[0.66, 1.0, 0.5],
													[0.0, 1.0, 0.5],
												],
								}}
								pathTracingSettings={roomConfig.PathTracer}
							/>
						)}

					{/* Per-particle vectors */}
					{roomConfig.VectorDisplay.vectors &&
						roomConfig.VectorDisplay.vectors.length > 0 &&
						perParticleVectors.length > 0 && (
							<PerParticleVectors
								vectors={perParticleVectors}
								arrowsConfig={{
									...roomConfig.VectorDisplay,
									colormap: vectorColormap,
								}}
								pathTracingSettings={roomConfig.PathTracer}
							/>
						)}

					{/* Simulation cell */}
					{currentFrame.cell.length > 0 &&
						roomConfig.Visualization.simulation_box && (
							<SimulationCell frame={currentFrame} colorMode={colorMode} />
						)}

					{/* Main particle visualization */}
					<ParticleInstances
						frame={currentFrame}
						setFrame={setCurrentFrame}
						selectedIds={selectedIds}
						setSelectedIds={setSelectedIds}
						isDrawing={isDrawing}
						setPoints={setPoints}
						setHoveredId={setHoveredId}
						sceneSettings={roomConfig.Particle}
						token={token}
						visibleIndices={undefined}
						highlight=""
						pathTracingSettings={roomConfig.PathTracer}
						roomLock={roomLock}
						editMode={editMode}
						setEditMode={setEditMode}
					/>

					{/* Hovered particle highlight */}
					{hoveredId !== -1 && (
						<ParticleInstances
							frame={currentFrame}
							setFrame={setCurrentFrame}
							selectedIds={selectedIds}
							setSelectedIds={setSelectedIds}
							isDrawing={isDrawing}
							setPoints={setPoints}
							setHoveredId={setHoveredId}
							sceneSettings={roomConfig.Particle}
							token={token}
							visibleIndices={hoveredId}
							highlight="backside"
							pathTracingSettings={roomConfig.PathTracer}
							roomLock={roomLock}
							editMode={editMode}
							setEditMode={setEditMode}
						/>
					)}

					{/* Selected particles highlight */}
					{selectedIds.size > 0 && (
						<>
							<ParticleInstances
								frame={currentFrame}
								setFrame={setCurrentFrame}
								selectedIds={selectedIds}
								setSelectedIds={setSelectedIds}
								isDrawing={isDrawing}
								setPoints={setPoints}
								setHoveredId={setHoveredId}
								sceneSettings={roomConfig.Particle}
								token={token}
								visibleIndices={selectedIds}
								highlight="selection"
								pathTracingSettings={roomConfig.PathTracer}
								roomLock={roomLock}
								editMode={editMode}
								setEditMode={setEditMode}
							/>
							<BondInstances
								frame={currentFrame}
								visibleIndices={selectedIds}
								highlight="selection"
								sceneSettings={roomConfig.Particle}
								pathTracingSettings={roomConfig.PathTracer}
							/>
						</>
					)}

					{/* Constrained particles highlight */}
					{constrainedAtomIds.size > 0 && (
						<ParticleInstances
							frame={currentFrame}
							setFrame={setCurrentFrame}
							selectedIds={selectedIds}
							setSelectedIds={setSelectedIds}
							isDrawing={isDrawing}
							setPoints={setPoints}
							setHoveredId={setHoveredId}
							sceneSettings={roomConfig.Particle}
							token={token}
							visibleIndices={constrainedAtomIds}
							highlight="constraint"
							pathTracingSettings={roomConfig.PathTracer}
							roomLock={roomLock}
							editMode={editMode}
							setEditMode={setEditMode}
						/>
					)}

					{/* Bonds */}
					<BondInstances
						frame={currentFrame}
						visibleIndices={undefined}
						sceneSettings={roomConfig.Particle}
						highlight=""
						pathTracingSettings={roomConfig.PathTracer}
					/>

					{/* Transform controls */}
					{selectedPoint && (
						<ControlsBuilder
							points={points}
							setPoints={setPoints}
							selectedPoint={selectedPoint}
							setSelectedPoint={setSelectedPoint}
						/>
					)}

					{/* 3D Lines */}
					{points.length > 0 && (
						<Line3D
							points={points}
							setPoints={setPoints}
							setSelectedPoint={setSelectedPoint}
							isDrawing={isDrawing}
							colorMode={colorMode}
							hoveredId={hoveredId}
							setIsDrawing={setIsDrawing}
							setLineLength={(length: number) => {
								// Line length is now calculated automatically via hook
								// This callback can be used for additional side effects if needed
							}}
						/>
					)}

					{/* Geometries */}
					<Geometries
						geometries={geometries}
						isDrawing={isDrawing}
						setPoints={setPoints}
						setHoveredId={setHoveredId}
					/>

					{/* Virtual Canvas for drawing */}
					<VirtualCanvas
						isDrawing={isDrawing}
						setPoints={setPoints}
						points={points}
						hoveredId={hoveredId}
						setHoveredId={setHoveredId}
					/>
				</Pathtracer>

				{/* Player component */}
				<Player
					playing={playing}
					step={step}
					setStep={setStep}
					fps={roomConfig.Camera.fps || 30}
					length={length}
					loop={roomConfig.Visualization.animation_loop || false}
					togglePlaying={setPlaying}
					selectedFrames={selectedFrames}
					isFrameRendering={isFrameRendering}
					setIsFrameRendering={setIsFrameRendering}
					frameRate={frameRate}
				/>
			</Canvas>

			{/* Virtual Canvas for 2D overlays - moved inside Canvas */}
		</div>
	);
};
