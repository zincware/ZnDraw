import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { getFrames } from "../../myapi/client";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect, useCallback } from "react";
import { useGeometryEditing } from "../../hooks/useGeometryEditing";
import { useGeometryPersistence } from "../../hooks/useGeometryPersistence";
import { renderMaterial } from "./materials";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import {
	processPositionAttribute,
	processRotationAttribute,
	processColorData,
	processScaleAttribute,
	getInstanceCount,
	validateArrayLengths,
	expandSharedColor,
	SELECTION_SCALE,
	HOVER_SCALE,
	type ScaleProp,
} from "../../utils/geometryData";
import {
	_vec3,
	_vec3_2,
	_euler,
	_matrix,
	_matrix2,
	_quat,
	_quat2,
	_color,
} from "../../utils/threeObjectPools";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";
import { useFrameKeys } from "../../hooks/useSchemas";

interface InteractionSettings {
	enabled: boolean;
	color: string;
	opacity: number;
}

interface ShapeData {
	vertices: [number, number][];
	position: string | number[][];
	rotation: string | number[][];
	color: string | string[];
	material: string;
	scale: ScaleProp;
	opacity: number;
	selecting: InteractionSettings;
	hovering: InteractionSettings;
	active?: boolean;
}

/**
 * Build a THREE.Shape from 2D vertices.
 * Connects vertices with straight lines and auto-closes.
 */
function buildShapeFromVertices(vertices: [number, number][]): THREE.Shape {
	const shape = new THREE.Shape();

	if (vertices.length === 0) return shape;

	// Start at first vertex
	shape.moveTo(vertices[0][0], vertices[0][1]);

	// Connect all subsequent vertices with lines
	for (let i = 1; i < vertices.length; i++) {
		shape.lineTo(vertices[i][0], vertices[i][1]);
	}

	// Auto-close the shape
	shape.closePath();

	return shape;
}

/**
 * Shape geometry component - renders 2D polygons in 3D space.
 * Supports instancing: one shape template rendered at multiple positions.
 *
 * Note: Vertex editing is not supported because Shape vertices are 2D
 * but transform controls operate in 3D.
 */
export default function Shape({
	data,
	geometryKey,
	pathtracingEnabled = false,
}: {
	data: ShapeData;
	geometryKey: string;
	pathtracingEnabled?: boolean;
}) {
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);

	// Merge with defaults from Pydantic (single source of truth)
	const fullData = getGeometryWithDefaults<ShapeData>(
		data,
		"Shape",
		geometryDefaults,
	);

	const {
		vertices,
		position: positionProp,
		rotation: rotationProp,
		color: colorProp,
		material,
		scale,
		opacity,
		selecting,
		hovering,
	} = fullData;

	const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
	const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
	const hoverMeshRef = useRef<THREE.Mesh | null>(null);

	const [instanceCount, setInstanceCount] = useState(0);

	// Use individual selectors to prevent unnecessary re-renders
	const currentFrame = useAppStore((state) => state.currentFrame);
	const frameCount = useAppStore((state) => state.frameCount);
	const roomId = useAppStore((state) => state.roomId);
	const userName = useAppStore((state) => state.userName);
	const selections = useAppStore((state) => state.selections);
	const updateSelections = useAppStore((state) => state.updateSelections);
	const hoveredGeometryInstance = useAppStore(
		(state) => state.hoveredGeometryInstance,
	);
	const setHoveredGeometryInstance = useAppStore(
		(state) => state.setHoveredGeometryInstance,
	);
	const setGeometryFetching = useAppStore((state) => state.setGeometryFetching);
	const removeGeometryFetching = useAppStore(
		(state) => state.removeGeometryFetching,
	);

	// Fetch frame keys to check if required data is available
	const { data: frameKeysData, isLoading: isLoadingKeys } = useFrameKeys(
		roomId!,
		currentFrame,
	);

	// Check if required keys are available
	const requiredKeys = useMemo(() => {
		const keys: string[] = [];
		if (typeof positionProp === "string") keys.push(positionProp);
		return keys;
	}, [positionProp]);

	const hasRequiredKeys = useMemo(() => {
		if (isLoadingKeys) return true;
		if (!frameKeysData?.keys) return true;
		const availableKeys = new Set(frameKeysData.keys);
		return requiredKeys.every((key) => availableKeys.has(key));
	}, [frameKeysData, requiredKeys, isLoadingKeys]);

	// Use geometry-specific selection
	const shapeSelection = selections[geometryKey] || [];
	const selectionSet = useMemo(() => new Set(shapeSelection), [shapeSelection]);
	const selectedIndices = useMemo(
		() => Array.from(selectionSet),
		[selectionSet],
	);
	const validSelectedIndices = useMemo(
		() => selectedIndices.filter((id) => id < instanceCount),
		[selectedIndices, instanceCount],
	);

	// Build THREE.Shape geometry from vertices
	const shapeGeometry = useMemo(() => {
		if (!vertices || vertices.length < 3) {
			return new THREE.ShapeGeometry(new THREE.Shape());
		}
		const shape = buildShapeFromVertices(vertices);
		return new THREE.ShapeGeometry(shape);
	}, [vertices]);

	// Individual queries for each attribute
	const { data: positionData, isFetching: isPositionFetching } = useQuery({
		queryKey: ["frame", roomId, currentFrame, positionProp],
		queryFn: ({ signal }: { signal: AbortSignal }) =>
			getFrames(roomId!, currentFrame, [positionProp as string], signal),
		enabled:
			!!roomId &&
			!!userName &&
			frameCount > 0 &&
			typeof positionProp === "string",
		placeholderData: keepPreviousData,
		retry: false,
	});

	const { data: colorData, isFetching: isColorFetching } = useQuery({
		queryKey: ["frame", roomId, currentFrame, colorProp],
		queryFn: ({ signal }: { signal: AbortSignal }) =>
			getFrames(roomId!, currentFrame, [colorProp as string], signal),
		enabled:
			!!roomId &&
			!!userName &&
			frameCount > 0 &&
			typeof colorProp === "string" &&
			shouldFetchAsFrameData(colorProp as string),
		placeholderData: keepPreviousData,
		retry: false,
	});

	const { data: rotationData, isFetching: isRotationFetching } = useQuery({
		queryKey: ["frame", roomId, currentFrame, rotationProp],
		queryFn: ({ signal }: { signal: AbortSignal }) =>
			getFrames(roomId!, currentFrame, [rotationProp as string], signal),
		enabled:
			!!roomId &&
			!!userName &&
			frameCount > 0 &&
			typeof rotationProp === "string",
		placeholderData: keepPreviousData,
		retry: false,
	});

	const { data: scaleData, isFetching: isScaleFetching } = useQuery({
		queryKey: ["frame", roomId, currentFrame, scale],
		queryFn: ({ signal }: { signal: AbortSignal }) =>
			getFrames(roomId!, currentFrame, [scale as string], signal),
		enabled:
			!!roomId && !!userName && frameCount > 0 && typeof scale === "string",
		placeholderData: keepPreviousData,
		retry: false,
	});

	// Check if any enabled query is still fetching
	const isFetching =
		(typeof positionProp === "string" && isPositionFetching) ||
		(typeof colorProp === "string" &&
			shouldFetchAsFrameData(colorProp as string) &&
			isColorFetching) ||
		(typeof rotationProp === "string" && isRotationFetching) ||
		(typeof scale === "string" && isScaleFetching);

	// Report fetching state to global store
	useEffect(() => {
		setGeometryFetching(geometryKey, isFetching);
	}, [geometryKey, isFetching, setGeometryFetching]);

	// Clean up fetching state on unmount
	useEffect(() => {
		return () => {
			removeGeometryFetching(geometryKey);
		};
	}, [geometryKey, removeGeometryFetching]);

	// Handle geometry editing with transform controls (for shape instances)
	const finalPositionData =
		typeof positionProp === "string"
			? positionData?.[positionProp]
			: positionProp;
	const finalRotationData =
		typeof rotationProp === "string"
			? rotationData?.[rotationProp]
			: rotationProp;
	const finalScaleData = typeof scale === "string" ? scaleData?.[scale] : scale;
	const scaleValue = finalScaleData || 1.0;

	useGeometryEditing(
		geometryKey,
		finalPositionData,
		finalRotationData,
		scaleValue,
		selectedIndices,
		"Shape",
		fullData,
		instanceCount,
	);

	// Handle persistence of local geometry changes to server
	useGeometryPersistence(geometryKey, "Shape");

	// Consolidated data processing and mesh update
	useEffect(() => {
		// When frameCount is 0, explicitly clear shapes (e.g., after del vis[:])
		if (frameCount === 0) {
			if (instanceCount !== 0) setInstanceCount(0);
			return;
		}

		if (isFetching) {
			return;
		}

		try {
			// --- Data Processing Step ---
			const fetchedPosition =
				typeof positionProp === "string"
					? positionData?.[positionProp as string]
					: undefined;
			const finalCount = getInstanceCount(positionProp, fetchedPosition);

			if (finalCount === 0) {
				if (instanceCount !== 0) setInstanceCount(0);
				return;
			}

			// Process all attributes
			const finalPositions = processPositionAttribute(
				positionProp,
				fetchedPosition,
			);

			const fetchedColor =
				typeof colorProp === "string"
					? colorData?.[colorProp as string]
					: undefined;
			const colorHexArray = processColorData(
				colorProp,
				fetchedColor,
				finalCount,
			);

			const fetchedRotation =
				typeof rotationProp === "string"
					? rotationData?.[rotationProp as string]
					: undefined;
			const finalRotations = processRotationAttribute(
				rotationProp,
				fetchedRotation,
				finalCount,
			);

			const fetchedScale =
				typeof scale === "string" ? scaleData?.[scale] : undefined;
			const { values: finalScales } = processScaleAttribute(
				scale,
				fetchedScale,
				finalCount,
			);

			// Handle shared color
			const finalColorHex = expandSharedColor(colorHexArray, finalCount);

			// --- Validation Step ---
			const isDataValid =
				validateArrayLengths(
					{
						positions: finalPositions,
						rotations: finalRotations,
						scales: finalScales,
					},
					{
						positions: finalCount * 3,
						rotations: finalCount * 3,
						scales: finalCount * 3,
					},
				) && finalColorHex.length === finalCount;

			if (!isDataValid) {
				console.error("Shape data is invalid or has inconsistent lengths.");
				if (instanceCount !== 0) setInstanceCount(0);
				return;
			}

			// --- Mesh Resizing Step ---
			if (instanceCount !== finalCount) {
				setInstanceCount(finalCount);
				return;
			}

			// --- Main Mesh Instance Update ---
			const mainMesh = mainMeshRef.current;
			if (!mainMesh) return;

			for (let i = 0; i < finalCount; i++) {
				const i3 = i * 3;
				_vec3.set(
					finalPositions[i3],
					finalPositions[i3 + 1],
					finalPositions[i3 + 2],
				);
				_euler.set(
					finalRotations[i3],
					finalRotations[i3 + 1],
					finalRotations[i3 + 2],
				);
				_quat.setFromEuler(_euler);
				_vec3_2.set(finalScales[i3], finalScales[i3 + 1], finalScales[i3 + 2]);
				_matrix.compose(_vec3, _quat, _vec3_2);
				mainMesh.setMatrixAt(i, _matrix);

				_color.set(finalColorHex[i]);
				mainMesh.setColorAt(i, _color);
			}

			mainMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
			mainMesh.instanceMatrix.needsUpdate = true;
			if (mainMesh.instanceColor) mainMesh.instanceColor.needsUpdate = true;

			mainMesh.computeBoundingBox();
			mainMesh.computeBoundingSphere();

			// --- Selection Mesh Update ---
			if (selecting.enabled && selectionMeshRef.current) {
				const selectionMesh = selectionMeshRef.current;
				validSelectedIndices.forEach((id, index) => {
					if (id >= finalCount) return;
					const i3 = id * 3;
					_vec3.set(
						finalPositions[i3],
						finalPositions[i3 + 1],
						finalPositions[i3 + 2],
					);
					_euler.set(
						finalRotations[i3],
						finalRotations[i3 + 1],
						finalRotations[i3 + 2],
					);
					_quat.setFromEuler(_euler);
					_vec3_2.set(
						finalScales[i3] * SELECTION_SCALE,
						finalScales[i3 + 1] * SELECTION_SCALE,
						finalScales[i3 + 2] * SELECTION_SCALE,
					);
					_matrix.compose(_vec3, _quat, _vec3_2);
					selectionMesh.setMatrixAt(index, _matrix);
				});
				selectionMesh.instanceMatrix.needsUpdate = true;
				selectionMesh.computeBoundingBox();
				selectionMesh.computeBoundingSphere();
			}
		} catch (error) {
			console.error("Error processing Shape data:", error);
			if (instanceCount !== 0) setInstanceCount(0);
		}
	}, [
		data,
		frameCount, // Watch frameCount to clear shapes when it becomes 0
		isFetching,
		positionData,
		colorData,
		rotationData,
		scaleData, // Add scaleData
		positionProp,
		colorProp,
		rotationProp,
		scale, // Add scale
		instanceCount,
		validSelectedIndices,
		selecting,
		geometryKey,
	]);

	// Separate effect for hover mesh updates
	useEffect(() => {
		if (!hovering?.enabled || !hoverMeshRef.current || !mainMeshRef.current)
			return;
		if (instanceCount === 0) return;

		const hoverMesh = hoverMeshRef.current;
		const mainMesh = mainMeshRef.current;

		if (
			hoveredGeometryInstance?.geometryKey === geometryKey &&
			hoveredGeometryInstance?.instanceId !== null &&
			hoveredGeometryInstance.instanceId < instanceCount
		) {
			hoverMesh.visible = true;

			mainMesh.getMatrixAt(hoveredGeometryInstance.instanceId, _matrix2);
			_matrix2.decompose(_vec3, _quat2, _vec3_2);

			hoverMesh.position.copy(_vec3);
			hoverMesh.quaternion.copy(_quat2);
			hoverMesh.scale.copy(_vec3_2).multiplyScalar(HOVER_SCALE);
		} else {
			hoverMesh.visible = false;
		}
	}, [hoveredGeometryInstance, instanceCount, hovering, geometryKey]);

	// Event handlers
	const onClickHandler = useCallback(
		(event: any) => {
			if (event.detail !== 1 || event.instanceId === undefined) return;
			event.stopPropagation();
			updateSelections(geometryKey, event.instanceId, event.shiftKey);
		},
		[updateSelections, geometryKey],
	);

	const onPointerEnterHandler = useCallback(
		(event: any) => {
			if (event.instanceId === undefined) return;
			event.stopPropagation();
			setHoveredGeometryInstance(geometryKey, event.instanceId);
		},
		[setHoveredGeometryInstance, geometryKey],
	);

	const onPointerOutHandler = useCallback(() => {
		setHoveredGeometryInstance(null, null);
	}, [setHoveredGeometryInstance]);

	if (!userName || !roomId) return null;

	// Don't render if geometry is disabled OR if required keys are not available
	if (fullData.active === false || !hasRequiredKeys) {
		return null;
	}

	return (
		<group>
			{/* Main instanced mesh */}
			<instancedMesh
				key={instanceCount}
				ref={mainMeshRef}
				args={[undefined, undefined, instanceCount]}
				visible={!pathtracingEnabled}
				onClick={
					!pathtracingEnabled && selecting.enabled ? onClickHandler : undefined
				}
				onPointerEnter={
					!pathtracingEnabled && hovering?.enabled
						? onPointerEnterHandler
						: undefined
				}
				onPointerOut={
					!pathtracingEnabled && hovering?.enabled
						? onPointerOutHandler
						: undefined
				}
			>
				<primitive object={shapeGeometry} attach="geometry" />
				{renderMaterial(material, opacity, undefined, THREE.DoubleSide)}
			</instancedMesh>

			{/* Selection mesh */}
			{!pathtracingEnabled && selecting.enabled && (
				<instancedMesh
					key={`selection-${validSelectedIndices.length}`}
					ref={selectionMeshRef}
					args={[undefined, undefined, validSelectedIndices.length]}
				>
					<primitive object={shapeGeometry} attach="geometry" />
					<meshBasicMaterial
						side={THREE.DoubleSide}
						transparent
						opacity={selecting.opacity}
						color={selecting.color}
					/>
				</instancedMesh>
			)}

			{/* Hover mesh */}
			{!pathtracingEnabled && hovering?.enabled && (
				<mesh ref={hoverMeshRef} visible={false}>
					<primitive object={shapeGeometry} attach="geometry" />
					<meshBasicMaterial
						side={THREE.DoubleSide}
						transparent
						opacity={hovering.opacity}
						color={hovering.color}
					/>
				</mesh>
			)}
		</group>
	);
}
