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
}: PlayerProps) => {
	useFrame(({ clock }) => {
		const a = clock.getElapsedTime();
		if (a > 1 / fps) {
			if (playing) {
				if (step === length - 1) {
					if (!loop) {
						setPlaying(!playing);
					} else {
						if (selectedFrames.indices.size > 0 && selectedFrames.active) {
							setStep(Math.min(...selectedFrames.indices));
						} else {
							setStep(0);
						}
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
						setStep(step + 1);
					}
				}
			}
			// reset the clock
			clock.start();
		}
	});
	return null;
};

const ParticleBondMaterial = ({
	highlight,
	material,
	color,
}: {
	highlight: string;
	material: string;
	color?: string;
}) => {
	if (highlight) {
		switch (highlight) {
			case "backside":
				return (
					<meshBasicMaterial side={THREE.BackSide} transparent opacity={0.8} />
				);
			case "frontside":
				return (
					<meshBasicMaterial side={THREE.FrontSide} transparent opacity={0.5} />
				);
			case "selection":
				return (
					<meshBasicMaterial side={THREE.FrontSide} transparent opacity={0.5} />
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
}) => {
	const meshRef = useRef<THREE.InstancedMesh | null>(null);

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
	const { selection_color, material, particle_size } = sceneSettings;

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
			const scaleVector = new THREE.Vector3();

			const visibleArray = Array.from(actualVisibleIndices);
			visibleArray.forEach((atomIdx, i) => {
				const position = positions[atomIdx];
				if (!position) {
					// if position was removed, this can happen and we skip it until next update
					return;
				}

				let radius;
				if (highlight === "backside") {
					radius = radii[atomIdx] * 1.25;
				} else if (highlight === "frontside") {
					radius = radii[atomIdx] * 1.01;
				} else if (highlight === "selection") {
					radius = radii[atomIdx] * 1.01;
				} else {
					radius = radii[atomIdx];
				}
				// Set position and scale for each instance
				const matrix = new THREE.Matrix4()
					.setPosition(position)
					.scale(scaleVector.set(radius, radius, radius));
				meshRef.current.setMatrixAt(i, matrix);
				color.set(highlight ? selection_color : colors[atomIdx]);
				meshRef.current.setColorAt(i, color);
			});

			// Mark instance matrices and colors for update
			meshRef.current.instanceMatrix.needsUpdate = true;
			meshRef.current.instanceColor.needsUpdate = true;
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

	if (highlight === "") {
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
				<ParticleBondMaterial highlight={highlight} material={material} />
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

	const { material, selection_color, bond_size } = sceneSettings;

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

	if (highlight === "") {
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
			<ParticleBondMaterial highlight={highlight} material={material} />
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
	frame: Frame | undefined;
	property: string;
	colorMode: string;
	arrowsConfig: any;
	pathTracingSettings: any;
}

export const PerParticleVectors: React.FC<PerParticleVectorsProps> = ({
	frame,
	property,
	colorMode,
	arrowsConfig,
	pathTracingSettings,
}) => {
	const [vectors, setVectors] = useState<
		{ start: THREE.Vector3; end: THREE.Vector3 }[]
	>([]);

	const LineColor = colorMode === "light" ? "#454b66" : "#f5fdc6";
	const LineWidth = 2;
	const LineScale = 1;

	const [colorRange, setColorRange] = useState<[number, number]>(
		arrowsConfig.colorrange,
	);

	useEffect(() => {
		if (arrowsConfig.normalize) {
			const max = Math.max(
				...vectors.map((vector) => vector.start.distanceTo(vector.end)),
			);
			setColorRange([0, max]);
		} else {
			setColorRange(arrowsConfig.colorrange);
		}
	}, [vectors, arrowsConfig.normalize, arrowsConfig.colorrange]);

	useEffect(() => {
		if (!frame) {
			setVectors([]);
			return;
		}
		let frameData;
		if (property in frame.calc) {
			frameData = frame.calc[property];
		} else if (property in frame.arrays) {
			frameData = frame.arrays[property];
		} else {
			console.error(`Property ${property} not found in frame`);
			setVectors([]);
			return;
		}
		if (frameData.length !== frame.positions.length) {
			console.error(
				`Length of property ${property} does not match the number of particles`,
			);
			setVectors([]);
			return;
		}

		console.log(`Property ${property} found in frame.calc`);
		const calculatedVectors = frameData.map((vector, i) => {
			const start = frame.positions[i];
			const end = start
				.clone()
				.add(
					new THREE.Vector3(vector[0], vector[1], vector[2]).multiplyScalar(
						LineScale,
					),
				);
			return { start, end };
		});
		setVectors(calculatedVectors);
	}, [frame, property]);

	const startMap = useMemo(
		() => vectors.map((vec) => vec.start.toArray()),
		[vectors],
	);
	const endMap = useMemo(
		() => vectors.map((vec) => vec.end.toArray()),
		[vectors],
	);

	return (
		<>
			<Arrows
				start={startMap}
				end={endMap}
				scale_vector_thickness={arrowsConfig.scale_vector_thickness}
				colormap={arrowsConfig.colormap}
				colorrange={colorRange}
				opacity={arrowsConfig.opacity}
				rescale={arrowsConfig.rescale}
				pathTracingSettings={pathTracingSettings}
			/>
		</>
	);
};
