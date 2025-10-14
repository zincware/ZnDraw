import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { getFrames } from "../../myapi/client";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect, useCallback } from "react";
import { renderMaterial } from "./materials";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import {
  processPositionAttribute,
  processRotationAttribute,
  processColorAttribute,
  processSize3D,
  getInstanceCount,
  validateArrayLengths,
  SELECTION_SCALE,
  HOVER_SCALE,
} from "../../utils/geometryData";
import { _vec3, _euler, _matrix, _color } from "../../utils/threeObjectPools";
import { convertInstancedMeshToMerged, disposeMesh } from "../../utils/convertInstancedMesh";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

interface InteractionSettings {
  enabled: boolean;
  color: string;
  opacity: number;
}

interface BoxData {
  position: string | number[][] | number[];
  size: string | number[][] | number[];
  color: string | number[][] | number[];
  rotation: string | number[][] | number[];
  material: string;
  scale: number;
  opacity: number;
  selecting: InteractionSettings;
  hovering: InteractionSettings;
}

export default function Box({
  data,
  geometryKey,
  pathtracingEnabled = false
}: {
  data: BoxData;
  geometryKey: string;
  pathtracingEnabled?: boolean;
}) {
  const { geometryDefaults } = useAppStore();

  // Merge with defaults from Pydantic (single source of truth)
  const fullData = getGeometryWithDefaults<BoxData>(data, "Box", geometryDefaults);

  const {
    position: positionProp,
    size: sizeProp,
    color: colorProp,
    rotation: rotationProp,
    material,
    scale,
    selecting,
    hovering,
  } = fullData;

  const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const hoverMeshRef = useRef<THREE.Mesh | null>(null);
  const mergedMeshRef = useRef<THREE.Mesh | null>(null);
  const [instanceCount, setInstanceCount] = useState(0);

  const { currentFrame, roomId, clientId, selections, updateSelections, hoveredGeometryInstance, setHoveredGeometryInstance, setDrawingPointerPosition, isDrawing, setDrawingIsValid, setGeometryFetching, removeGeometryFetching, requestPathtracingUpdate } = useAppStore();

  // Use geometry-specific selection
  const boxSelection = selections[geometryKey] || [];
  const selectionSet = useMemo(() => new Set(boxSelection), [boxSelection]);
  const selectedIndices = useMemo(() => Array.from(selectionSet), [selectionSet]);
  const validSelectedIndices = useMemo(
    () => selectedIndices.filter((id) => id < instanceCount),
    [selectedIndices, instanceCount]
  );

  const boxScale = scale || 1.0;

  // Individual queries for each attribute
  const { data: positionData, isFetching: isPositionFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, positionProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [positionProp as string], signal),
    enabled: !!roomId && !!clientId && typeof positionProp === "string",
    placeholderData: keepPreviousData,
  });

  const { data: sizeData, isFetching: isSizeFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, sizeProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [sizeProp as string], signal),
    enabled: !!roomId && !!clientId && typeof sizeProp === "string",
    placeholderData: keepPreviousData,
  });

  const { data: colorData, isFetching: isColorFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, colorProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [colorProp as string], signal),
    enabled: !!roomId && !!clientId && typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string),
    placeholderData: keepPreviousData,
  });

  const { data: rotationData, isFetching: isRotationFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, rotationProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [rotationProp as string], signal),
    enabled: !!roomId && !!clientId && typeof rotationProp === "string",
    placeholderData: keepPreviousData,
  });

  // Check if any enabled query is still fetching
  const isFetching =
    (typeof positionProp === "string" && isPositionFetching) ||
    (typeof sizeProp === "string" && isSizeFetching) ||
    (typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string) && isColorFetching) ||
    (typeof rotationProp === "string" && isRotationFetching);

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

  // Consolidated data processing and mesh update
  useEffect(() => {
    if (isFetching) {
      return; // Wait for all enabled queries to complete
    }

    try {

      // --- Data Processing Step ---
      const fetchedPosition = typeof positionProp === 'string' ? positionData?.[positionProp as string] : undefined;
      const finalCount = getInstanceCount(positionProp, fetchedPosition);

      if (finalCount === 0) {
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      // Process all attributes using new specialized functions
      const finalPositions = processPositionAttribute(positionProp, fetchedPosition);

      const fetchedColor = typeof colorProp === 'string' ? colorData?.[colorProp as string] : undefined;
      const finalColors = processColorAttribute(colorProp, fetchedColor, finalCount);

      const fetchedSize = typeof sizeProp === 'string' ? sizeData?.[sizeProp as string] : undefined;
      const finalSizes = processSize3D(sizeProp, fetchedSize, finalCount);

      const fetchedRotation = typeof rotationProp === 'string' ? rotationData?.[rotationProp as string] : undefined;
      const finalRotations = processRotationAttribute(rotationProp, fetchedRotation, finalCount);

      // --- Validation Step ---
      const isDataValid = validateArrayLengths(
        { positions: finalPositions, colors: finalColors, sizes: finalSizes, rotations: finalRotations },
        { positions: finalCount * 3, colors: finalCount * 3, sizes: finalCount * 3, rotations: finalCount * 3 }
      );

      if (!isDataValid) {
        console.error("Box data is invalid or has inconsistent lengths.");
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
        _vec3.set(finalPositions[i3], finalPositions[i3 + 1], finalPositions[i3 + 2]);
        _euler.set(finalRotations[i3], finalRotations[i3 + 1], finalRotations[i3 + 2]);
        const width = finalSizes[i3] * boxScale;
        const height = finalSizes[i3 + 1] * boxScale;
        const depth = finalSizes[i3 + 2] * boxScale;
        _matrix.compose(_vec3, new THREE.Quaternion().setFromEuler(_euler), new THREE.Vector3(width, height, depth));
        mainMesh.setMatrixAt(i, _matrix);

        // Use original colors (no inline selection/hover colors)
        _color.setRGB(finalColors[i3], finalColors[i3 + 1], finalColors[i3 + 2]);
        mainMesh.setColorAt(i, _color);
      }

      mainMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
      mainMesh.instanceMatrix.needsUpdate = true;
      if (mainMesh.instanceColor) mainMesh.instanceColor.needsUpdate = true;

      // --- Selection Mesh Update ---
      if (selecting.enabled && selectionMeshRef.current) {
        const selectionMesh = selectionMeshRef.current;
        validSelectedIndices.forEach((id, index) => {
          if (id >= finalCount) return;
          const i3 = id * 3;
          _vec3.set(finalPositions[i3], finalPositions[i3 + 1], finalPositions[i3 + 2]);
          _euler.set(finalRotations[i3], finalRotations[i3 + 1], finalRotations[i3 + 2]);
          const width = finalSizes[i3] * boxScale * SELECTION_SCALE;
          const height = finalSizes[i3 + 1] * boxScale * SELECTION_SCALE;
          const depth = finalSizes[i3 + 2] * boxScale * SELECTION_SCALE;
          _matrix.compose(_vec3, new THREE.Quaternion().setFromEuler(_euler), new THREE.Vector3(width, height, depth));
          selectionMesh.setMatrixAt(index, _matrix);
        });
        selectionMesh.instanceMatrix.needsUpdate = true;
      }

    } catch (error) {
      console.error("Error processing Box data:", error);
      if (instanceCount !== 0) setInstanceCount(0);
    }
  }, [
    isFetching,
    positionData,
    sizeData,
    colorData,
    rotationData,
    positionProp,
    sizeProp,
    colorProp,
    rotationProp,
    instanceCount,
    boxScale,
    validSelectedIndices,
    selecting,
    geometryKey,
  ]);

  // Separate effect for hover mesh updates - doesn't trigger data reprocessing
  useEffect(() => {
    if (!hovering?.enabled || !hoverMeshRef.current || !mainMeshRef.current) return;
    if (instanceCount === 0) return;

    const hoverMesh = hoverMeshRef.current;
    const mainMesh = mainMeshRef.current;

    // Only show hover if it's for this geometry
    if (hoveredGeometryInstance?.geometryKey === geometryKey &&
        hoveredGeometryInstance?.instanceId !== null &&
        hoveredGeometryInstance.instanceId < instanceCount) {
      hoverMesh.visible = true;

      // Get transform from main mesh
      const matrix = new THREE.Matrix4();
      mainMesh.getMatrixAt(hoveredGeometryInstance.instanceId, matrix);

      const position = new THREE.Vector3();
      const quaternion = new THREE.Quaternion();
      const scale = new THREE.Vector3();
      matrix.decompose(position, quaternion, scale);

      // Apply hover scale
      hoverMesh.position.copy(position);
      hoverMesh.quaternion.copy(quaternion);
      hoverMesh.scale.copy(scale).multiplyScalar(HOVER_SCALE);
    } else {
      hoverMesh.visible = false;
    }
  }, [hoveredGeometryInstance, instanceCount, hovering, geometryKey]);

  // Convert instanced mesh to merged mesh for path tracing
  useEffect(() => {
    if (!pathtracingEnabled) {
      if (mergedMeshRef.current) {
        disposeMesh(mergedMeshRef.current);
        mergedMeshRef.current = null;
      }
      return;
    }

    if (!mainMeshRef.current || instanceCount === 0) return;

    if (mergedMeshRef.current) {
      disposeMesh(mergedMeshRef.current);
    }

    const mergedMesh = convertInstancedMeshToMerged(mainMeshRef.current);
    mergedMeshRef.current = mergedMesh;

    requestPathtracingUpdate();

    return () => {
      if (mergedMeshRef.current) {
        disposeMesh(mergedMeshRef.current);
        mergedMeshRef.current = null;
      }
    };
  }, [
    pathtracingEnabled,
    instanceCount,
    geometryKey,
    requestPathtracingUpdate,
  ]);

  // Shared geometry for all boxes
  const mainGeometry = useMemo(() => {
    return new THREE.BoxGeometry(1, 1, 1);
  }, []);

  const onClickHandler = useCallback((event: any) => {
    if (event.detail !== 1 || event.instanceId === undefined) return;
    event.stopPropagation();
    updateSelections(geometryKey, event.instanceId, event.shiftKey);
  }, [updateSelections, geometryKey]);

  const onPointerMoveHandler = useCallback((event: any) => {
    if (event.instanceId === undefined) return;
    event.stopPropagation();
    if (isDrawing) {
      setDrawingPointerPosition(event.point);
    }
  }, [isDrawing, setDrawingPointerPosition]);

  const onPointerEnterHandler = useCallback((event: any) => {
    if (event.instanceId === undefined) return;
    event.stopPropagation();
    setHoveredGeometryInstance(geometryKey, event.instanceId);
    setDrawingIsValid(true);
  }, [setHoveredGeometryInstance, setDrawingIsValid, geometryKey]);

  const onPointerOutHandler = useCallback(() => {
    setHoveredGeometryInstance(null, null);
    setDrawingIsValid(false);
  }, [setHoveredGeometryInstance, setDrawingIsValid]);

  if (!clientId || !roomId) return null;

  return (
    <group>
      {/* Main instanced mesh - visible when NOT pathtracing */}
      <instancedMesh
        key={instanceCount}
        ref={mainMeshRef}
        args={[undefined, undefined, instanceCount]}
        visible={!pathtracingEnabled}
        onClick={!pathtracingEnabled && selecting.enabled ? onClickHandler : undefined}
        onPointerEnter={!pathtracingEnabled && hovering?.enabled ? onPointerEnterHandler : undefined}
        onPointerMove={!pathtracingEnabled && hovering?.enabled ? onPointerMoveHandler : undefined}
        onPointerOut={!pathtracingEnabled && hovering?.enabled ? onPointerOutHandler : undefined}
      >
        <primitive object={mainGeometry} attach="geometry" />
        {renderMaterial(material, data.opacity)}
      </instancedMesh>

      {/* Selection mesh - only when NOT pathtracing */}
      {!pathtracingEnabled && selecting.enabled && (
        <instancedMesh
          key={`selection-${validSelectedIndices.length}`}
          ref={selectionMeshRef}
          args={[undefined, undefined, validSelectedIndices.length]}
        >
          <primitive object={mainGeometry} attach="geometry" />
          <meshBasicMaterial
            side={THREE.FrontSide}
            transparent
            opacity={selecting.opacity}
            color={selecting.color}
          />
        </instancedMesh>
      )}

      {/* Hover mesh - only when NOT pathtracing */}
      {!pathtracingEnabled && hovering?.enabled && (
        <mesh ref={hoverMeshRef} visible={false}>
          <primitive object={mainGeometry} attach="geometry" />
          <meshBasicMaterial
            side={THREE.BackSide}
            transparent
            opacity={hovering.opacity}
            color={hovering.color}
          />
        </mesh>
      )}

      {/* Merged mesh - visible when pathtracing */}
      {pathtracingEnabled && mergedMeshRef.current && (
        <primitive object={mergedMeshRef.current} />
      )}
    </group>
  );
}
