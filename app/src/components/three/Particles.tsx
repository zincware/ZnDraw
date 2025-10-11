import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { getFrames } from "../../myapi/client";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect, useCallback } from "react";
import { renderMaterial } from "./materials";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import {
  type DataProp,
  processNumericAttribute,
  processColorAttribute,
  getInstanceCount,
  validateArrayLengths,
  SELECTION_SCALE,
  HOVER_SCALE,
} from "../../utils/geometryData";
import { _vec3, _matrix, _color } from "../../utils/threeObjectPools";

interface InteractionSettings {
  enabled: boolean;
  color: string;
  opacity: number;
}

interface SphereData {
  position: string | number[][] | number[];
  color: string | number[][] | number[];
  radius: string | number[] | number;
  material: string;
  resolution: number;
  scale: number;
  opacity: number;
  selecting: InteractionSettings;
  hovering: InteractionSettings;
}

// Reusable THREE objects imported from threeObjectPools

export default function Sphere({ data, geometryKey }: { data: SphereData; geometryKey: string }) {
  const {
    position: positionProp,
    color: colorProp,
    radius: radiusProp,
    material,
    resolution,
    scale,
    selecting,
    hovering,
    opacity,
  } = data;

  const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const hoverMeshRef = useRef<THREE.Mesh | null>(null);
  const [hoveredInstanceId, setHoveredInstanceId] = useState<number | null>(null);
  const [instanceCount, setInstanceCount] = useState(0);

  const { currentFrame, roomId, clientId, selection, updateSelection, setDrawingPointerPosition, isDrawing, setDrawingIsValid, setGeometryFetching, removeGeometryFetching } = useAppStore();

  const selectionSet = useMemo(() => new Set(selection || []), [selection]);
  const selectedIndices = useMemo(() => Array.from(selectionSet), [selectionSet]);
  const validSelectedIndices = useMemo(
    () => selectedIndices.filter((id) => id < instanceCount),
    [selectedIndices, instanceCount]
  );

  const particleResolution = resolution || 8;
  const particleScale = scale || 1.0;
  // console.log("queryKey", ["frame", roomId, currentFrame, positionProp]);
  // Individual queries for each attribute - enables perfect cross-component caching
  const { data: positionData, isFetching: isPositionFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, positionProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [positionProp as string], signal),
    enabled: !!roomId && !!clientId && typeof positionProp === "string",
    placeholderData: keepPreviousData,
  });

  const { data: colorData, isFetching: isColorFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, colorProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [colorProp as string], signal),
    enabled: !!roomId && !!clientId && typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string),
    placeholderData: keepPreviousData,
  });

  const { data: radiusData, isFetching: isRadiusFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, radiusProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [radiusProp as string], signal),
    enabled: !!roomId && !!clientId && typeof radiusProp === "string",
    placeholderData: keepPreviousData,
  });

  // Check if any enabled query is still fetching
  const isFetching =
    (typeof positionProp === "string" && isPositionFetching) ||
    (typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string) && isColorFetching) ||
    (typeof radiusProp === "string" && isRadiusFetching);

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

      // Process all attributes
      const finalPositions = processNumericAttribute(positionProp, fetchedPosition, finalCount);

      const fetchedColor = typeof colorProp === 'string' ? colorData?.[colorProp as string] : undefined;
      const finalColors = processColorAttribute(colorProp, fetchedColor, finalCount);

      const fetchedRadius = typeof radiusProp === 'string' ? radiusData?.[radiusProp as string] : undefined;
      const finalRadii = processNumericAttribute(radiusProp, fetchedRadius, finalCount);

      // --- Validation Step ---
      const isDataValid = validateArrayLengths(
        { positions: finalPositions, colors: finalColors, radii: finalRadii },
        { positions: finalCount * 3, colors: finalCount * 3, radii: finalCount }
      );

      if (!isDataValid) {
        console.error("Sphere/Particles data is invalid or has inconsistent lengths.");
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
        const r = finalRadii[i] * particleScale;
        _matrix.identity().setPosition(_vec3).scale(_vec3.set(r, r, r));
        mainMesh.setMatrixAt(i, _matrix);
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
          const r = finalRadii[id] * particleScale * SELECTION_SCALE;
          _matrix.identity().setPosition(_vec3).scale(_vec3.set(r, r, r));
          selectionMesh.setMatrixAt(index, _matrix);
        });
        selectionMesh.instanceMatrix.needsUpdate = true;
      }

      // --- Hover Mesh Update ---
      if (hovering?.enabled && hoverMeshRef.current) {
        const hoverMesh = hoverMeshRef.current;
        if (hoveredInstanceId !== null && hoveredInstanceId < finalCount) {
          hoverMesh.visible = true;
          const i3 = hoveredInstanceId * 3;
          _vec3.set(finalPositions[i3], finalPositions[i3 + 1], finalPositions[i3 + 2]);
          const r = finalRadii[hoveredInstanceId] * particleScale * HOVER_SCALE;
          hoverMesh.position.copy(_vec3);
          hoverMesh.scale.set(r, r, r);
        } else {
          hoverMesh.visible = false;
        }
      }
    } catch (error) {
      console.error("Error processing Sphere/Particles data:", error);
      if (instanceCount !== 0) setInstanceCount(0);
    }
  }, [
    isFetching,
    positionData,
    colorData,
    radiusData,
    positionProp,
    colorProp,
    radiusProp,
    instanceCount,
    particleScale,
    validSelectedIndices,
    selecting,
    hovering,
    hoveredInstanceId,
  ]);


  const mainGeometry = useMemo(() => {
    return new THREE.SphereGeometry(1, particleResolution, particleResolution);
  }, [particleResolution]);


  const onClickHandler = useCallback((event: any) => {
    if (event.detail !== 1 || event.instanceId === undefined) return;
    event.stopPropagation();
    updateSelection(event.instanceId, event.shiftKey);
  }, [updateSelection]);

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
    setHoveredInstanceId(event.instanceId);
    setDrawingIsValid(true);
  }, [setDrawingIsValid]);

  const onPointerOutHandler = useCallback(() => { setHoveredInstanceId(null); setDrawingIsValid(false); }, [setDrawingIsValid]);

  if (!clientId || !roomId) return null;

  return (
    <group>
      {/* Main mesh */}
      <instancedMesh
        key={instanceCount}
        ref={mainMeshRef}
        args={[undefined, undefined, instanceCount]}
        onClick={selecting.enabled ? onClickHandler : undefined}
        onPointerEnter={hovering?.enabled ? onPointerEnterHandler : undefined}
        onPointerMove={hovering?.enabled ? onPointerMoveHandler : undefined}
        onPointerOut={hovering?.enabled ? onPointerOutHandler : undefined}
      >
        <primitive object={mainGeometry} attach="geometry" />
        {renderMaterial(material, data.opacity)}
      </instancedMesh>

      {/* Selection mesh */}
      {selecting.enabled && (
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

      {/* Hover mesh */}
      {hovering?.enabled && (
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
    </group>
  );
}