import React, { useMemo } from "react";
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
import { Line3D } from "../lines";
import ControlsBuilder from "../transforms";
import VectorField from "../vectorfield";
import { getCentroid } from "../particlesEditor";
import * as THREE from "three";

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
	} = useAppContext();

	// Calculate derived values
	const centroid = useMemo(() => {
		if (selectedIds && selectedIds.size > 0) {
			return getCentroid(currentFrame.positions, selectedIds);
		}
		return getCentroid(currentFrame.positions, new Set());
	}, [selectedIds, currentFrame.positions]);

	const lineLength = useMemo(() => {
		if (points.length >= 2) {
			const start = points[points.length - 2];
			const end = points[points.length - 1];
			return start.distanceTo(end);
		}
		return 0;
	}, [points]);

	return (
		<div
			className="canvas-container edit-mode-border"
			onDragOver={onDragOver}
			onDrop={onDrop}
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
							<Environment preset={roomConfig.PathTracer.environment} />
						)}

					{/* Lighting */}
					{roomConfig.Visualization.directional_light &&
						(roomConfig.Visualization.floor ? (
							<>
								<Floor colorMode={colorMode} roomConfig={roomConfig} />
								<directionalLight
									position={[0, 100, 0]}
									intensity={1.0}
									castShadow
									shadow-mapSize-width={roomConfig.Camera.camera_far * 10}
									shadow-mapSize-height={roomConfig.Camera.camera_far * 10}
									shadow-camera-near={10}
									shadow-camera-far={800}
									shadow-camera-left={-1 * roomConfig.Camera.camera_far}
									shadow-camera-right={roomConfig.Camera.camera_far}
									shadow-camera-top={roomConfig.Camera.camera_far}
									shadow-camera-bottom={-1 * roomConfig.Camera.camera_far}
								/>
							</>
						) : (
							<directionalLight position={[0, 100, 0]} intensity={1.0} />
						))}

					{/* Vector field */}
					{roomConfig.VectorDisplay.vectorfield &&
						currentFrame.vectors !== undefined && (
							<VectorField
								vectors={currentFrame.vectors}
								arrowsConfig={roomConfig.Arrows}
								pathTracingSettings={roomConfig.PathTracer}
							/>
						)}

					{/* Per-particle vectors */}
					{roomConfig.VectorDisplay.property &&
						roomConfig.VectorDisplay.property !== "none" && (
							<PerParticleVectors
								frame={currentFrame}
								step={step}
								property={roomConfig.VectorDisplay.property}
								colorMode={colorMode}
								arrowsConfig={roomConfig.VectorDisplay}
								pathTracingSettings={roomConfig.PathTracer}
								token={token}
							/>
						)}

					{/* Simulation cell */}
					{currentFrame.cell.length > 0 && (
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
							setLineLength={() => {}} // TODO: Add lineLength management
						/>
					)}

					{/* Geometries */}
					<Geometries
						geometries={geometries}
						isDrawing={isDrawing}
						setPoints={setPoints}
						setHoveredId={setHoveredId}
					/>
				</Pathtracer>

				{/* Player component */}
				<Player
					playing={playing}
					step={step}
					setStep={setStep}
					fps={roomConfig.Visualization.fps}
					length={length}
					loop={roomConfig.Visualization.loop}
					togglePlaying={setPlaying}
					selectedFrames={selectedFrames}
				/>
			</Canvas>

			{/* Virtual Canvas for 2D overlays - moved inside Canvas */}
		</div>
	);
};
