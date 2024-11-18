import { Pathtracer } from "@react-three/gpu-pathtracer";
import React, { useState, useEffect, useRef } from "react";
import {
	setupBookmarks,
	setupCamera,
	setupConfig,
	setupFigures,
	setupFrames,
	setupGeometries,
	setupMessages,
	setupPoints,
	setupSelection,
	setupStep,
} from "./components/api";
import { Geometries } from "./components/geometries";
import HeadBar from "./components/headbar";
import {
	BondInstances,
	ParticleInstances,
	PerParticleVectors,
	SimulationCell,
} from "./components/particles";
import { type Frame, Frames, Player } from "./components/particles";
import { getCentroid } from "./components/particlesEditor";
import FrameProgressBar from "./components/progressbar";
import Sidebar from "./components/sidebar";
import { client, socket } from "./socket";
import "./App.css";
import * as znsocket from "znsocket";
import { Plotting } from "./components/plotting";

import {
	Box,
	CameraControls,
	Environment,
	OrbitControls,
	OrthographicCamera,
	PerspectiveCamera,
	Sphere,
	TrackballControls,
	type TransformControls,
} from "@react-three/drei";
import { Canvas, useFrame } from "@react-three/fiber";
import { Button, Form, InputGroup } from "react-bootstrap";
import * as THREE from "three";
import CameraAndControls from "./components/cameraAndControls";
import { Floor } from "./components/floor";
import { Line3D, VirtualCanvas } from "./components/lines";
import { ParticleInfoOverlay, SceneInfoOverlay } from "./components/overlays";
import ControlsBuilder from "./components/transforms";
import { useColorMode } from "./components/utils";
import type { IndicesState } from "./components/utils";
import VectorField from "./components/vectorfield";

export default function App() {
	// const [isConnected, setIsConnected] = useState(socket.connected);
	// const [fooEvents, setFooEvents] = useState([]);

	const [connected, setConnected] = useState<boolean>(false);
	const [currentFrame, setCurrentFrame] = useState<Frame>({
		arrays: { colors: [], radii: [] },
		calc: null,
		cell: [],
		connectivity: [],
		info: null,
		numbers: [],
		pbc: [],
		positions: [],
		vectors: [],
	});
	const [playing, setPlaying] = useState<boolean>(false);
	const [length, setLength] = useState<number>(0);
	// updated via sockets
	const [step, setStep] = useState<number>(0);
	const [selectedFrames, setSelectedFrames] = useState<IndicesState>({
		active: true,
		indices: new Set<number>(),
	});
	const [selectedIds, setSelectedIds] = useState<Set<number>>(new Set());
	const [bookmarks, setBookmarks] = useState<any>({}); // {name: [step, ...]
	const [points, setPoints] = useState<THREE.Vector3[]>([]);

	const [isDrawing, setIsDrawing] = useState<boolean>(false);
	const [selectedPoint, setSelectedPoint] = useState<THREE.Vector3 | null>(
		null,
	);
	const [roomName, setRoomName] = useState<string>("");
	const [geometries, setGeometries] = useState<any>([]);
	const [cameraAndControls, setCameraAndControls] = useState<any>({
		camera: new THREE.Vector3(0, 0, 0),
		target: new THREE.Vector3(0, 0, 0),
	});
	// TODO: initial values are wrong for orbitcontrolstarget and camperaPosition
	// todo give to particles and bonds
	const [colorMode, handleColorMode] = useColorMode();
	const [hoveredId, setHoveredId] = useState<number>(-1);
	// UPDATE THESE using `vis.config` on the Python side
	const [roomConfig, setRoomConfig] = useState({
		Arrows: {
			colormap: [
				[-0.5, 0.9, 0.5],
				[0.0, 0.9, 0.5],
			],
			normalize: true,
			colorrange: [0, 1.0],
			scale_vector_thickness: false,
			opacity: 1.0,
		},
		Particle: {
			particle_size: 1.0,
			bond_size: 1.0,
			material: "MeshStandardMaterial",
			selection_color: "#ffa500",
		},
		Visualization: {
			simulation_box: true,
			floor: false,
			frame_update: true,
			animation_loop: false,
		},
		Camera: {
			camera: "PerspectiveCamera",
			camera_near: 0.1,
			camera_far: 300,
			crosshair: false,
			synchronize_camera: true,
			fps: 30,
			controls: "OrbitControls",
		},
		PathTracer: {
			enabled: false,
			environment: "studio",
			metalness: 0.7,
			roughness: 0.2,
			clearcoat: 0.0,
			clearcoatRoughness: 0.0,
		},
		VectorDisplay: {
			vectorfield: true,
			vectors: "",
			vector_scale: 1.0,
		},
	});

	const [isAuthenticated, setIsAuthenticated] = useState<boolean>(true);
	const [roomLock, setRoomLock] = useState<boolean>(false);

	// QUEUES
	// TODO: fix
	const [modifierQueue, setModifierQueue] = useState<number>(-1);

	// extension UI elements
	const [tutorialURL, setTutorialURL] = useState<string>("");
	const [showSiMGen, setShowSiMGen] = useState<boolean>(false);

	const [cursorPosition, setCursorPosition] = useState({ x: 0, y: 0 });
	const [lineLength, setLineLength] = useState<number>(0);
	const [showParticleInfo, setShowParticleInfo] = useState<boolean>(false);
	const [addPlotsWindow, setAddPlotsWindow] = useState<number>(0); // make this bool!
	const [updatedPlotsList, setUpdatedPlotsList] = useState<string[]>([]);
	const [messages, setMessages] = useState<string[]>([]);

	const [token, setToken] = useState<string>("");
	setupConfig(token, setRoomConfig);
	setupBookmarks(token, setBookmarks, bookmarks);
	setupPoints(token, setPoints, points);
	setupSelection(token, setSelectedIds, selectedIds);
	setupStep(token, setStep, step);
	setupCamera(
		token,
		cameraAndControls,
		setCameraAndControls,
		roomConfig.Camera.synchronize_camera,
	);
	setupFrames(
		token,
		step,
		setCurrentFrame,
		currentFrame,
		setLength,
		setStep,
		roomConfig.Visualization.frame_update,
	);
	setupFigures(token, setUpdatedPlotsList);
	setupGeometries(token, setGeometries, geometries);
	setupMessages(token, setMessages, messages);

	useEffect(() => {
		function onConnect() {
			setConnected(true);
			socket.emit(
				"webclient:connect",
				(data: { name: string; room: string; authenticated: boolean }) => {
					setRoomName(data.room);
					setIsAuthenticated(data.authenticated);
				},
			);
			console.log("connected");

			// get lock state
			socket.emit("room:lock:get", (data: boolean) => {
				setRoomLock(data);
			});

			socket.emit("room:token:get", (data: string) => {
				setToken(data);
			});
		}

		function onDisconnect() {
			setConnected(false);
		}

		function onTutorialURL(data: string) {
			setTutorialURL(data);
		}
		function onShowSiMGen(data: boolean) {
			setShowSiMGen(data);
		}

		function onRoomLockSet(locked: boolean) {
			setRoomLock(locked);
		}

		socket.on("connect", onConnect);
		socket.on("disconnect", onDisconnect);
		socket.on("tutorial:url", onTutorialURL);
		socket.on("showSiMGen", onShowSiMGen);
		socket.on("room:lock:set", onRoomLockSet);

		return () => {
			socket.off("connect", onConnect);
			socket.off("disconnect", onDisconnect);
			socket.off("tutorial:url", onTutorialURL);
			socket.off("showSiMGen", onShowSiMGen);
			socket.off("room:lock:set", onRoomLockSet);
		};
	}, []);

	useEffect(() => {
		// page initialization

		const handleKeyDown = (event: KeyboardEvent) => {
			// if canvas is not focused, don't do anything
			if (document.activeElement !== document.body) {
				return;
			}
			if (event.key === "ArrowRight") {
				setPlaying(false);
				if (event.shiftKey) {
					// Jump to the bookmark after the current step
					const bookmarkKeys = Object.keys(bookmarks)
						.map(Number)
						.sort((a, b) => a - b);
					const nextBookmark = bookmarkKeys.find((key) => key > step);

					if (nextBookmark !== undefined) {
						setStep(nextBookmark);
					}
				} else {
					if (selectedFrames.indices.size > 0 && selectedFrames.active) {
						const nextFrame = Array.from(selectedFrames.indices).find(
							(frame) => frame > step,
						);
						if (nextFrame) {
							setStep(nextFrame);
						} else {
							setStep(Math.min(...selectedFrames.indices));
						}
					} else {
						setStep((prevStep) => (prevStep + 1 < length ? prevStep + 1 : 0));
					}
				}
			} else if (event.key === "ArrowLeft") {
				setPlaying(false);
				if (event.shiftKey) {
					// Jump to the bookmark before the current step
					const bookmarkKeys = Object.keys(bookmarks)
						.map(Number)
						.sort((a, b) => b - a);
					const previousBookmark = bookmarkKeys.find((key) => key < step);

					if (previousBookmark !== undefined) {
						setStep(previousBookmark);
					}
				} else {
					// Move to the previous step, or wrap around to the end
					// check if selectedFrames length is greater than 0, then only jump
					// between selectedFrames
					if (selectedFrames.indices.size > 0 && selectedFrames.active) {
						const previousFrame = Array.from(selectedFrames.indices)
							.reverse()
							.find((frame) => frame < step);
						if (previousFrame) {
							setStep(previousFrame);
						} else {
							setStep(Math.max(...selectedFrames.indices));
						}
					} else {
						setStep((prevStep) =>
							prevStep - 1 >= 0 ? prevStep - 1 : length - 1,
						);
					}
				}
			} else if (event.key === "ArrowUp") {
				// jump 10 percent, or to the end
				setPlaying(false);
				const newStep = Math.min(step + Math.floor(length / 10), length - 1);
				setStep(newStep);
			} else if (event.key === "ArrowDown") {
				// jump 10 percent, or to the beginning
				setPlaying(false);
				const newStep = Math.max(step - Math.floor(length / 10), 0);
				setStep(newStep);
			} else if (event.key === " ") {
				// backspace
				// updateLength();
				setPlaying((prev) => !prev);
				if (step === length - 1) {
					setStep(0);
				}
			} else if (event.key === "x") {
				setIsDrawing((prev) => !prev);
			} else if (event.key === "i") {
				setShowParticleInfo((prev) => !prev);
			} else if (event.key === "b") {
				setBookmarks((prev) => {
					const newBookmarks = { ...prev };
					newBookmarks[step] = `Frame ${step}`;
					return newBookmarks;
				});
			} else if (event.key === "a") {
				if (event.ctrlKey) {
					setSelectedIds(
						new Set([
							...Array.from(
								{ length: currentFrame.positions.length },
								(_, i) => i,
							),
						]),
					);
				}
			} else if (event.key === "Backspace" || event.key === "Delete") {
				// check if shift is pressed
				if (event.shiftKey) {
					if (selectedPoint !== null) {
						const newPoints = points.filter(
							(point) => point.distanceTo(selectedPoint) > 0.1,
						);
						setSelectedPoint(null);
						setPoints(newPoints);
					} else if (points.length > 0) {
						// pop last point from points
						setSelectedPoint(null);
						setPoints(points.slice(0, points.length - 1));
					}
				} else {
					if (selectedIds.size > 0) {
						const queue = new znsocket.Dict({
							client: client,
							key: `queue:${token}:modifier`,
						});
						queue.Delete = {};
						socket.emit("room:worker:run");
					}
				}
			}
		};

		// Add the event listener
		window.addEventListener("keydown", handleKeyDown);

		// Clean up the event listener on unmount
		return () => {
			window.removeEventListener("keydown", handleKeyDown);
		};
	}, [
		length,
		step,
		points,
		selectedPoint,
		bookmarks,
		currentFrame,
		selectedIds,
	]);

	// reduce selection, if selected points is reduced
	useEffect(() => {
		if (selectedIds.size > 0) {
			const newSelectedIds = new Set(
				Array.from(selectedIds).filter(
					(id) => id < currentFrame.positions.length,
				),
			);
			// if the selection is reduced, update the selection
			if (newSelectedIds.size < selectedIds.size) {
				setSelectedIds(newSelectedIds);
			}
		}
	}, [currentFrame, selectedIds]);

	useEffect(() => {
		const updateCursorPosition = (event) => {
			setCursorPosition({ x: event.clientX, y: event.clientY });
		};

		window.addEventListener("mousemove", updateCursorPosition);

		return () => {
			window.removeEventListener("mousemove", updateCursorPosition);
		};
	}, []);

	const onDragOver = (event) => {
		event.preventDefault();
	};

	const onDrop = async (event) => {
		event.preventDefault();

		const file = event.dataTransfer.files[0];
		if (!file) {
			console.error("No file was dropped");
			return;
		}

		try {
			const arrayBuffer = await file.arrayBuffer();
			const content = new Uint8Array(arrayBuffer);

			// send the file to the server
			socket.emit("room:upload:file", {
				content: Array.from(content),
				filename: file.name,
			});
		} catch (error) {
			console.error("Error reading file:", error);
		}
	};

	const onPointerMissed = () => {
		setSelectedPoint(null);
		setSelectedIds(new Set());
	};

	return (
		<>
			<div className="canvas-container" onDragOver={onDragOver} onDrop={onDrop}>
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

						{roomConfig.Visualization.floor ? (
							<>
								<Floor colorMode={colorMode} roomConfig={roomConfig} />
								<directionalLight
									position={[0, 100, 0]}
									intensity={1.0}
									castShadow
									shadow-mapSize-width={roomConfig.Camera.camera_far * 10} // Adjust the width of the shadow map
									shadow-mapSize-height={roomConfig.Camera.camera_far * 10} // Adjust the height of the shadow map
									shadow-camera-near={10} // Adjust the near clipping plane of the shadow camera
									shadow-camera-far={800} // Adjust the far clipping plane of the shadow camera
									shadow-camera-left={-1 * roomConfig.Camera.camera_far} // Set the left boundary for the shadow camera frustum
									shadow-camera-right={roomConfig.Camera.camera_far} // Set the right boundary for the shadow camera frustum
									shadow-camera-top={roomConfig.Camera.camera_far} // Set the top boundary for the shadow camera frustum
									shadow-camera-bottom={-1 * roomConfig.Camera.camera_far} // Set the bottom boundary for the shadow camera frustum
								/>
							</>
						) : (
							<directionalLight position={[0, 100, 0]} intensity={1.0} />
						)}

						{roomConfig.VectorDisplay.vectorfield &&
							currentFrame.vectors !== undefined && (
								<VectorField
									vectors={currentFrame.vectors}
									pathTracingSettings={roomConfig.PathTracer}
									arrowsConfig={{
										rescale: roomConfig.VectorDisplay.vector_scale,
										...roomConfig.Arrows,
									}}
								/>
							)}
						<ParticleInstances
							frame={currentFrame}
							selectedIds={selectedIds}
							setSelectedIds={setSelectedIds}
							isDrawing={isDrawing}
							setPoints={setPoints}
							setHoveredId={setHoveredId}
							sceneSettings={roomConfig.Particle}
							token={token}
							highlight=""
							visibleIndices={undefined}
							setFrame={setCurrentFrame}
							pathTracingSettings={roomConfig.PathTracer}
						/>
						{!roomConfig.PathTracer.enabled && (
							<>
								<ParticleInstances
									frame={currentFrame}
									selectedIds={selectedIds}
									setSelectedIds={setSelectedIds}
									isDrawing={isDrawing}
									setPoints={setPoints}
									setHoveredId={setHoveredId}
									sceneSettings={roomConfig.Particle}
									token={token}
									visibleIndices={hoveredId}
									highlight={"backside"}
									setFrame={setCurrentFrame}
								/>
								<ParticleInstances
									frame={currentFrame}
									selectedIds={selectedIds}
									setSelectedIds={setSelectedIds}
									isDrawing={isDrawing}
									setPoints={setPoints}
									setHoveredId={setHoveredId}
									sceneSettings={roomConfig.Particle}
									token={token}
									visibleIndices={selectedIds}
									highlight={"selection"}
									setFrame={setCurrentFrame}
								/>
								<ParticleInstances
									frame={currentFrame}
									selectedIds={selectedIds}
									setSelectedIds={setSelectedIds}
									isDrawing={isDrawing}
									setPoints={setPoints}
									setHoveredId={setHoveredId}
									sceneSettings={roomConfig.Particle}
									token={token}
									visibleIndices={
										new Set(currentFrame.constraints?.[0]?.indices)
									}
									highlight={"constraint"}
									setFrame={setCurrentFrame}
								/>
								<BondInstances
									frame={currentFrame}
									visibleIndices={selectedIds}
									highlight="selection"
									sceneSettings={roomConfig.Particle}
								/>
							</>
						)}
						<BondInstances
							frame={currentFrame}
							visibleIndices={undefined}
							highlight=""
							sceneSettings={roomConfig.Particle}
							pathTracingSettings={roomConfig.PathTracer}
						/>
						{roomConfig.Visualization.simulation_box &&
							!roomConfig.PathTracer.enabled && (
								<SimulationCell frame={currentFrame} colorMode={colorMode} />
							)}
						<Player
							playing={playing}
							togglePlaying={setPlaying}
							step={step}
							setStep={setStep}
							fps={roomConfig.Camera.fps}
							loop={roomConfig.Visualization.animation_loop}
							length={length}
							selectedFrames={selectedFrames}
						/>
						<Line3D
							points={points}
							setPoints={setPoints}
							setSelectedPoint={setSelectedPoint}
							isDrawing={isDrawing}
							colorMode={colorMode}
							hoveredId={hoveredId}
							setIsDrawing={setIsDrawing}
							setLineLength={setLineLength}
						/>
						<ControlsBuilder
							points={points}
							setPoints={setPoints}
							selectedPoint={selectedPoint}
							setSelectedPoint={setSelectedPoint}
						/>
						<Geometries
							geometries={geometries}
							isDrawing={isDrawing}
							setHoveredId={setHoveredId}
							setPoints={setPoints}
						/>
						<VirtualCanvas
							setPoints={setPoints}
							isDrawing={isDrawing}
							points={points}
							hoveredId={hoveredId}
							setHoveredId={setHoveredId}
						/>
						{roomConfig.VectorDisplay.vectors[0] &&
							roomConfig.VectorDisplay.vectors.map((vector) => (
								<PerParticleVectors
									frame={currentFrame}
									property={vector}
									colorMode={colorMode}
									arrowsConfig={{
										rescale: roomConfig.VectorDisplay.vector_scale,
										...roomConfig.Arrows,
									}}
									pathTracingSettings={roomConfig.PathTracer}
									key={vector}
								/>
							))}
					</Pathtracer>
				</Canvas>
			</div>
			<div className="App">
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
					length={length}
					step={step}
					setStep={setStep}
					bookmarks={bookmarks}
					setBookmarks={setBookmarks}
					selectedFrames={selectedFrames}
					setSelectedFrames={setSelectedFrames}
					connected={connected}
				/>
				<Plotting
					setStep={setStep}
					setSelectedFrames={setSelectedFrames}
					addPlotsWindow={addPlotsWindow}
					setSelectedIds={setSelectedIds}
					step={step}
					updatedPlotsList={updatedPlotsList}
					token={token}
				/>
				{showParticleInfo && (
					<>
						<ParticleInfoOverlay
							show={hoveredId !== null || isDrawing}
							info={{
								...(hoveredId !== null && {
									"Particle ID": hoveredId,
									"Atomic Number": currentFrame.numbers[hoveredId],
								}),
								...(isDrawing && { Line: `${lineLength.toFixed(2)} Å` }),
							}}
							position={cursorPosition}
						/>
						<SceneInfoOverlay
							frame={currentFrame}
							setShowParticleInfo={setShowParticleInfo}
						/>
					</>
				)}
			</div>
		</>
	);
}
