import type React from "react";
import { useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import * as znsocket from "znsocket";
import { client, socket } from "../socket";

import { Line } from "@react-three/drei";
import { useFrame } from "@react-three/fiber";
import Arrows from "./meshes";
import type { IndicesState } from "./utils";

import { ParticleControls } from "./particlesEditor";
import { useMergedMesh } from "./utils/mergeInstancedMesh";

export interface Frame {
	arrays: { colors: Array<string>; radii: Array<number> };
	calc: any;
	cell: number[][];
	connectivity: Array<[number, number, number]>;
	info: any;
	numbers: number[];
	pbc: boolean[];
	positions: THREE.Vector3[]; // only number[][] | before being mapped immediately
	vectors: [number, number, number][][];
	constraints: any[];
}

export interface Frames {
	[key: number]: { _type: string; value: Frame };
}

interface PlayerProps {
	playing: boolean;
	step: number;
	setStep: (step: number) => void;
	fps: number;
	length: number;
	loop: boolean;
	togglePlaying: (playing: boolean) => void;
	selectedFrames: IndicesState;
	isFrameRendering: boolean;
	setIsFrameRendering: (rendering: boolean) => void;
	frameRate: number;
}

export const Player = ({
	playing,
	step,
	setStep,
	fps,
	length,
	loop,
	togglePlaying: setPlaying,
	selectedFrames,
	isFrameRendering,
	setIsFrameRendering,
	frameRate,
}: PlayerProps) => {
	const lastUpdateTime = useRef(0);
	const frameTime = 1 / fps; // 30 FPS hardcoded for now

	// Handle play state changes - if starting to play on the last frame, jump to first frame
	useEffect(() => {
		if (playing) {
			if (selectedFrames.indices.size > 0 && selectedFrames.active) {
				// Selection mode - check if we're at the last selected frame
				const selectedFramesList = Array.from(selectedFrames.indices).sort((a, b) => a - b);
				const lastSelectedFrame = selectedFramesList[selectedFramesList.length - 1];
				if (step >= lastSelectedFrame) {
					setStep(selectedFramesList[0]); // Jump to first selected frame
				}
			} else {
				// Normal mode - check if we're at the last frame
				if (step >= length - 1) {
					setStep(0);
				}
			}
		}
	}, [playing, step, length, setStep, selectedFrames]);

	useFrame(({ clock }) => {
		if (playing && !isFrameRendering) {
			// check if the difference is greater than frameTime
			if (clock.getElapsedTime() - lastUpdateTime.current > frameTime) {
				lastUpdateTime.current = clock.getElapsedTime();

				// Mark frame as rendering
				setIsFrameRendering(true);

				setStep((prevStep) => {
					let nextStep = prevStep;

					// Check if frame selection is active
					if (selectedFrames.indices.size > 0 && selectedFrames.active) {
						// Playing within selected frames only
						const selectedFramesList = Array.from(selectedFrames.indices).sort((a, b) => a - b);
						const currentIndex = selectedFramesList.indexOf(prevStep);
						
						if (currentIndex !== -1) {
							// Current frame is in selection, move to next selected frame
							let nextIndex = currentIndex;
							
							// Apply frame rate (skip selected frames)
							for (let i = 0; i < frameRate && nextIndex < selectedFramesList.length - 1; i++) {
								nextIndex++;
							}
							
							if (nextIndex < selectedFramesList.length) {
								nextStep = selectedFramesList[nextIndex];
							} else {
								// Reached end of selected frames
								if (loop) {
									// Loop back to first selected frame
									nextStep = selectedFramesList[0];
								} else {
									// Stop playing at last selected frame
									nextStep = selectedFramesList[selectedFramesList.length - 1];
									setPlaying(false);
								}
							}
						} else {
							// Current frame is not in selection, find next selected frame
							const nextSelectedFrame = selectedFramesList.find(frame => frame > prevStep);
							if (nextSelectedFrame !== undefined) {
								nextStep = nextSelectedFrame;
							} else if (loop) {
								// Wrap to first selected frame
								nextStep = selectedFramesList[0];
							} else {
								// No more selected frames, stop playing
								setPlaying(false);
								nextStep = prevStep;
							}
						}
					} else {
						// Normal playback through all frames
						if (prevStep < length - 1) {
							nextStep = Math.min(prevStep + frameRate, length - 1);
							// Check if we've reached the last frame
							if (nextStep >= length - 1) {
								if (loop) {
									// If looping is enabled, continue playing (will wrap on next iteration)
									// Don't stop playing here
								} else {
									// If looping is disabled, stop playing when we reach the last frame
									setPlaying(false);
								}
							}
						} else if (prevStep >= length - 1) {
							// We're at or past the last frame
							if (loop) {
								// If looping is enabled, wrap to the first frame and continue playing
								nextStep = 0;
							} else {
								// If looping is disabled, stay at the last frame and stop playing
								nextStep = prevStep;
								setPlaying(false);
							}
						}
					}

					// Frame will be marked as finished after rendering completes
					return nextStep;
				});
			}
		}
	});

	return null;
};

const ParticleBondMaterial = ({
	highlight,
	material,
	color,
	hover_opacity,
	selection_opacity,
}: {
	highlight: string;
	material: string;
	color?: string;
	hover_opacity?: number;
	selection_opacity?: number;
}) => {
	if (highlight) {
		switch (highlight) {
			case "backside":
				return (
					<meshBasicMaterial
						side={THREE.BackSide}
						transparent
						opacity={hover_opacity}
					/>
				);
			case "selection":
				return (
					<meshBasicMaterial
						side={THREE.FrontSide}
						transparent
						opacity={selection_opacity}
					/>
				);
			case "constraint":
				return (
					<meshBasicMaterial side={THREE.FrontSide} color="red" wireframe />
				);
			default:
				return null;
		}
	}

	switch (material) {
		case "MeshPhysicalMaterial":
			return (
				<meshPhysicalMaterial
					roughness={0.3}
					reflectivity={0.4}
					side={THREE.FrontSide}
					color={color}
				/>
			);
		case "MeshToonMaterial":
			return <meshToonMaterial side={THREE.FrontSide} color={color} />;
		case "MeshStandardMaterial":
			return <meshStandardMaterial side={THREE.FrontSide} color={color} />;
		case "MeshBasicMaterial":
			return <meshBasicMaterial side={THREE.FrontSide} color={color} />;
		default:
			return null;
	}
};

export const ParticleInstances = ({
	frame,
	setFrame,
	selectedIds,
	setSelectedIds,
	isDrawing,
	setPoints,
	setHoveredId,
	sceneSettings,
	token,
	visibleIndices = undefined,
	highlight = "",
	pathTracingSettings,
	roomLock,
	editMode = "",
	setEditMode = () => {},
}: {
	frame: Frame;
	setFrame: (frame: Frame) => void;
	selectedIds: Set<number>;
	setSelectedIds: any;
	isDrawing: boolean;
	setPoints: any;
	setHoveredId: any;
	sceneSettings: any;
	token: string;
	visibleIndices: Set<number> | undefined | number;
	highlight: string;
	pathTracingSettings: any;
	roomLock: boolean;
	editMode?: string; // Make editMode optional
	setEditMode?: (mode: string) => void; // Make setEditMode optional
}) => {
	const meshRef = useRef<THREE.InstancedMesh | null>(null);


	useEffect(() => {
		console.log("Updating frame:", frame);
	}, [frame]);
	const actualVisibleIndices = useMemo(() => {
		if (typeof visibleIndices === "number") {
			if (visibleIndices === -1) {
				// -1 means no hover
				return new Set();
			}
			return new Set([visibleIndices]);
		}
		return (
			visibleIndices ??
			new Set(Array.from({ length: frame.numbers.length }, (_, i) => i))
		);
	}, [visibleIndices, frame.numbers.length]);

	const { colors, radii } = frame.arrays;
	const positions = frame.positions;
	const {
		selection_color,
		material,
		particle_size,
		hover_opacity,
		selection_opacity,
	} = sceneSettings;

	const geometry = useMemo(() => {
		const _geometry = new THREE.SphereGeometry(1, 32, 32);
		_geometry.scale(particle_size, particle_size, particle_size);
		if (pathTracingSettings?.enabled) {
			_geometry.scale(0, 0, 0);
		}
		return _geometry;
	}, [pathTracingSettings?.enabled, particle_size]);

	const instancedGeometry = useMemo(() => {
		return new THREE.SphereGeometry(1, 32, 32);
	}, []);

	useEffect(() => {
		if (meshRef.current && actualVisibleIndices.size > 0) {
			const color = new THREE.Color(selection_color);
			const matrix = new THREE.Matrix4(); // Reuse matrix object

			const visibleArray = Array.from(actualVisibleIndices);
			const highlightColor = highlight ? selection_color : null;
			
			// Pre-calculate scale multipliers to avoid repeated conditionals
			const scaleMultiplier = highlight === "backside" ? 1.25 : 
								   highlight === "selection" ? 1.01 : 1.0;

			for (let i = 0; i < visibleArray.length; i++) {
				const atomIdx = visibleArray[i];
				const position = positions[atomIdx];
				if (!position) {
					// if position was removed, this can happen and we skip it until next update
					continue;
				}

				const radius = radii[atomIdx] * scaleMultiplier;
				
				// Optimize matrix operations by setting elements directly
				matrix.makeScale(radius, radius, radius);
				matrix.setPosition(position);
				meshRef.current.setMatrixAt(i, matrix);
				
				// Set color efficiently
				color.set(highlightColor || colors[atomIdx] || "#ffffff");
				meshRef.current.setColorAt(i, color);
			}

			// Mark instance matrices and colors for update
			meshRef.current.instanceMatrix.needsUpdate = true;
			if (meshRef.current.instanceColor) {
				meshRef.current.instanceColor.needsUpdate = true;
			}
		}
	}, [
		positions,
		colors,
		radii,
		actualVisibleIndices,
		selection_color,
		highlight,
		geometry,
	]);

	const handlePointerOver = (event) => {
		if (highlight !== "") {
			return;
		}
		event.stopPropagation();
		setHoveredId(event.instanceId);
		// detect shift and control key being pressed at the same time
		if (event.shiftKey && event.ctrlKey) {
			selectedIds.add(event.instanceId);
			setSelectedIds(new Set(selectedIds));
		}
	};

	const handlePointerOut = (event) => {
		if (highlight !== "") {
			return;
		}
		event.stopPropagation();
		setHoveredId(-1);
	};

	const handleClicked = (event) => {
		if (event.detail !== 1) {
			return; // only handle single clicks
		}
		if (highlight !== "") {
			return;
		}

		event.stopPropagation();
		if (!event.shiftKey) {
			if (selectedIds.has(event.instanceId)) {
				setSelectedIds(new Set());
			} else {
				setSelectedIds(new Set([event.instanceId]));
			}
		} else {
			if (selectedIds.has(event.instanceId)) {
				selectedIds.delete(event.instanceId);
			} else {
				selectedIds.add(event.instanceId);
			}
			setSelectedIds(new Set(selectedIds));
		}
	};

	const handleDoubleClick = (event) => {
		const queue = new znsocket.Dict({
			client: client,
			key: `queue:${token}:selection`,
		});
		queue.ConnectedParticles = {};
		socket.emit("room:worker:run");
		event.stopPropagation();
	};

	const handlePointerMove = (event) => {
		if (!isDrawing || highlight) return;
		event.stopPropagation();

		const hoverPoint = event.intersections.find((i) => i.object.visible)?.point;
		if (hoverPoint) {
			setPoints((prevPoints: THREE.Vector3[]) => [
				...prevPoints.slice(0, -1),
				hoverPoint,
			]);
		}
	};

	if (highlight === "" && pathTracingSettings?.enabled) {
		useMergedMesh(meshRef, instancedGeometry, pathTracingSettings, [
			frame,
			visibleIndices,
		]);
	}

	return (
		<>
			{highlight === "selection" && (
				<ParticleControls
					frame={frame}
					selectedIds={selectedIds}
					setFrame={setFrame}
					roomLock={roomLock}
					editMode={editMode}
					setEditMode={setEditMode}
				/>
			)}
			<instancedMesh
				ref={meshRef}
				args={[geometry, null, actualVisibleIndices.size]}
				onPointerOver={handlePointerOver}
				onPointerOut={handlePointerOut}
				onPointerMove={highlight === "" ? handlePointerMove : undefined}
				onClick={handleClicked}
				onDoubleClick={handleDoubleClick}
				castShadow
				frustumCulled={false}
			>
				<ParticleBondMaterial
					highlight={highlight}
					material={material}
					hover_opacity={hover_opacity}
					selection_opacity={selection_opacity}
				/>
			</instancedMesh>
		</>
	);
};

export const BondInstances = ({
	frame,
	visibleIndices = undefined,
	highlight = "",
	sceneSettings,
	pathTracingSettings,
}: {
	frame: Frame;
	sceneSettings: any;
	visibleIndices: Set<number> | undefined;
	highlight: string;
	pathTracingSettings: any;
}) => {
	const meshRef = useRef<THREE.InstancedMesh | null>(null);

	const {
		material,
		selection_color,
		bond_size,
		hover_opacity,
		selection_opacity,
	} = sceneSettings;

	const actualVisibleConnectivity = useMemo(() => {
		if (!visibleIndices) {
			return frame.connectivity;
		}
		// find the subset of frame.connectivity where one of the particles are in visibleIndices
		return frame.connectivity.filter(([a, b]) => {
			return visibleIndices?.has(a) && visibleIndices?.has(b);
		});
	}, [visibleIndices, frame.connectivity]);

	const geometry = useMemo(() => {
		const _geometry = new THREE.CylinderGeometry(0.14, 0.14, 1, 32, 1, true);
		_geometry.translate(0, 0.5, 0);
		_geometry.rotateX(Math.PI / 2);

		if (!pathTracingSettings?.enabled) {
			_geometry.scale(bond_size, bond_size, 0.5);
		} else {
			// set to zero to make invisible
			_geometry.scale(0, 0, 0.5);
		}
		return _geometry;
	}, [pathTracingSettings?.enabled, bond_size]);

	const instancedGeometry = useMemo(() => {
		const _geometry = new THREE.CylinderGeometry(0.14, 0.14, 1, 32, 1, true);
		_geometry.translate(0, 0.5, 0);
		_geometry.rotateX(Math.PI / 2);

		_geometry.scale(1, 1, 0.5);
		return _geometry;
	}, []);

	useEffect(() => {
		if (meshRef.current && actualVisibleConnectivity.length > 0) {
			const color = new THREE.Color();
			const matrix = new THREE.Matrix4();
			const up = new THREE.Vector3(0, 1, 0);
			const direction = new THREE.Vector3();
			const scale = new THREE.Vector3(1, 1, 1);

			const createTransformationMatrix = (
				posA: THREE.Vector3,
				posB: THREE.Vector3,
			) => {
				direction.subVectors(posB, posA).normalize();
				const distance = posA.distanceTo(posB);
				scale.set(1, 1, distance);
				if (highlight === "selection") {
					scale.multiplyScalar(1.01);
				}

				matrix.lookAt(posA, posB, up);
				matrix.scale(scale);
				matrix.setPosition(posA.clone().lerp(posB, 0.5));

				return matrix.clone(); // Clone to avoid overwriting
			};

			actualVisibleConnectivity.forEach(([a, b], i) => {
				const posA = frame.positions[a] as THREE.Vector3;
				const posB = frame.positions[b] as THREE.Vector3;
				if (!posA || !posB) {
					// console.error("Connected particles not found");
					return;
				}

				// Set matrix and color for the bond from A to B
				meshRef.current?.setMatrixAt(
					i * 2,
					createTransformationMatrix(posA, posB),
				);
				if (highlight === "selection") {
					color.set(selection_color);
				} else {
					color.setStyle(frame.arrays.colors[a]);
				}
				meshRef.current?.setColorAt(i * 2, color);

				// Set matrix and color for the bond from B to A
				meshRef.current?.setMatrixAt(
					i * 2 + 1,
					createTransformationMatrix(posB, posA),
				);
				if (highlight === "selection") {
					color.set(selection_color);
				} else {
					color.setStyle(frame.arrays.colors[b]);
				}
				meshRef.current?.setColorAt(i * 2 + 1, color);
			});

			meshRef.current.instanceMatrix.needsUpdate = true;
			meshRef.current.instanceColor.needsUpdate = true;
		}
	}, [frame, actualVisibleConnectivity, selection_color, geometry]);

	if (highlight === "" && pathTracingSettings?.enabled) {
		useMergedMesh(meshRef, instancedGeometry, pathTracingSettings, [
			frame,
			visibleIndices,
		]);
	}

	return (
		<instancedMesh
			ref={meshRef}
			args={[geometry, undefined, actualVisibleConnectivity.length * 2]}
			castShadow
			// receiveShadow
		>
			<ParticleBondMaterial
				highlight={highlight}
				material={material}
				hover_opacity={hover_opacity}
				selection_opacity={selection_opacity}
			/>
		</instancedMesh>
	);
};

export const SimulationCell = ({
	frame,
	colorMode,
}: {
	frame: Frame;
	colorMode: string;
}) => {
	const [lineColor, setLineColor] = useState("black");

	useEffect(() => {
		if (colorMode === "light") {
			setLineColor("black");
		} else {
			setLineColor("white");
		}
	}, [colorMode]);

	const vertices = useMemo(() => {
		const cell = frame.cell;
		if (cell.length !== 3) {
			console.error("Invalid cell dimensions");
			return;
		}

		const origin = new THREE.Vector3(0, 0, 0);
		const vectors = cell.map((row) => new THREE.Vector3(...row));

		// Create the vertices of the box
		const v = [
			origin,
			vectors[0],
			vectors[1],
			vectors[0].clone().add(vectors[1]),
			vectors[2],
			vectors[0].clone().add(vectors[2]),
			vectors[1].clone().add(vectors[2]),
			vectors[0].clone().add(vectors[1]).add(vectors[2]),
		];

		return [
			[
				// Bottom face
				v[0],
				v[1],
				v[1],
				v[3],
				v[3],
				v[2],
				v[2],
				v[0],
			],
			[
				// Top face
				v[4],
				v[5],
				v[5],
				v[7],
				v[7],
				v[6],
				v[6],
				v[4],
			],
			// Vertical lines

			[v[0], v[4]],
			[v[1], v[5]],
			[v[2], v[6]],
			[v[3], v[7]],
		];
	}, [frame.cell]);

	return (
		<>
			{vertices && (
				<>
					<Line points={vertices[0]} color={lineColor} lineWidth={2} />
					<Line points={vertices[1]} color={lineColor} lineWidth={2} />
					<Line points={vertices[2]} color={lineColor} lineWidth={2} />
					<Line points={vertices[3]} color={lineColor} lineWidth={2} />
					<Line points={vertices[4]} color={lineColor} lineWidth={2} />
					<Line points={vertices[5]} color={lineColor} lineWidth={2} />
				</>
			)}
		</>
	);
};

interface PerParticleVectorsProps {
	vectors: { start: THREE.Vector3; end: THREE.Vector3; vectorType: string }[];
	arrowsConfig: any;
	pathTracingSettings: any;
}

export const PerParticleVectors: React.FC<PerParticleVectorsProps> = ({
	vectors,
	arrowsConfig,
	pathTracingSettings,
}) => {
	if (vectors.length === 0) {
		return null;
	}

	// Group vectors by type
	const vectorsByType = useMemo(() => {
		const grouped: Record<string, { start: THREE.Vector3; end: THREE.Vector3 }[]> = {};
		for (const vector of vectors) {
			if (!grouped[vector.vectorType]) {
				grouped[vector.vectorType] = [];
			}
			grouped[vector.vectorType].push({
				start: vector.start,
				end: vector.end,
			});
		}
		return grouped;
	}, [vectors]);

	// Helper function to convert hex color to HSL colormap
	const hexToHSLColormap = (hexColor: string): [number, number, number][] => {
		// Simple hex to HSL conversion for single color
		const hex = hexColor.replace('#', '');
		const r = parseInt(hex.substr(0, 2), 16) / 255;
		const g = parseInt(hex.substr(2, 2), 16) / 255;
		const b = parseInt(hex.substr(4, 2), 16) / 255;
		
		const max = Math.max(r, g, b);
		const min = Math.min(r, g, b);
		const diff = max - min;
		const l = (max + min) / 2;
		
		let h = 0;
		let s = 0;
		
		if (diff !== 0) {
			s = l > 0.5 ? diff / (2 - max - min) : diff / (max + min);
			
			switch (max) {
				case r: h = (g - b) / diff + (g < b ? 6 : 0); break;
				case g: h = (b - r) / diff + 2; break;
				case b: h = (r - g) / diff + 4; break;
			}
			h /= 6;
		}
		
		return [[h, s, l]];
	};

	return (
		<>
			{Object.entries(vectorsByType).map(([vectorType, typeVectors]) => {
				const vectorColor = arrowsConfig.vector_colors?.[vectorType] || "#ff0000";
				const colormap = hexToHSLColormap(vectorColor);
				
				const startMap = typeVectors.map((vec) => vec.start.toArray());
				const endMap = typeVectors.map((vec) => vec.end.toArray());
				
				// Calculate color range for this vector type
				const distances = typeVectors.map((vector) => vector.start.distanceTo(vector.end));
				const colorRange: [number, number] = arrowsConfig.normalize 
					? [0, Math.max(...distances)]
					: (Array.isArray(arrowsConfig.colorrange) ? arrowsConfig.colorrange : [0, 1]);

				return (
					<Arrows
						key={vectorType}
						start={startMap}
						end={endMap}
						scale_vector_thickness={arrowsConfig.scale_vector_thickness || false}
						colormap={colormap}
						colorrange={colorRange}
						opacity={arrowsConfig.opacity || 1.0}
						rescale={arrowsConfig.rescale}
						pathTracingSettings={pathTracingSettings}
					/>
				);
			})}
		</>
	);
};
