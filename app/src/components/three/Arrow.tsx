import * as THREE from "three";
import { useAppStore } from "../../store";
import { useRef, useMemo, useEffect, useState, useCallback } from "react";
import { BufferGeometryUtils } from "three/examples/jsm/Addons.js";
import { renderMaterial } from "./materials";
import { getFrames } from "../../myapi/client";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import {
  type DataProp,
  processNumericAttribute,
  processColorAttribute,
  getInstanceCount,
  validateArrayLengths,
  SELECTION_COLOR,
} from "../../utils/geometryData";
import { _vec3, _vec3_2, _vec3_3, _quat, _matrix, _color } from "../../utils/threeObjectPools";
import { convertInstancedMeshToMerged, disposeMesh } from "../../utils/convertInstancedMesh";

interface ArrowProps {
  position: DataProp;
  direction: DataProp;
  color: DataProp;
  radius: DataProp;
  scale: DataProp;
  material: string;
  geometryKey: string;
}

// Arrow-specific reusable THREE objects
const _arrowUp = new THREE.Vector3(0, 1, 0);

/**
 * Creates a standard arrow geometry with a total height of 1,
 * with its base at the origin, pointing along the positive Y-axis.
 */
function createArrowMesh() {
  const cylinderRadius = 0.04;
  const cylinderHeight = 0.6;
  const coneRadius = 0.1;
  const coneHeight = 0.4;

  const cylinderGeometry = new THREE.CylinderGeometry(
    cylinderRadius,
    cylinderRadius,
    cylinderHeight,
    32
  );
  const coneGeometry = new THREE.ConeGeometry(coneRadius, coneHeight, 32);

  // Position geometries so the base is at (0,0,0) and it extends up to a total height of 1.0
  cylinderGeometry.translate(0, cylinderHeight / 2, 0);
  coneGeometry.translate(0, cylinderHeight + coneHeight / 2, 0);

  const arrowGeometry = BufferGeometryUtils.mergeGeometries([
    cylinderGeometry,
    coneGeometry,
  ]);
  return arrowGeometry;
}

export default function Arrow({
  data,
  geometryKey,
  pathtracingEnabled = false
}: {
  data: ArrowProps;
  geometryKey: string;
  pathtracingEnabled?: boolean;
}) {
  const { position, direction, color, radius, scale, material } = data;
  const instancedMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const mergedMeshRef = useRef<THREE.Mesh | null>(null);
  const [instanceCount, setInstanceCount] = useState(0);

  const { currentFrame, roomId, clientId, selections, updateSelections, setGeometryFetching, removeGeometryFetching, requestPathtracingUpdate } = useAppStore();

  // Use geometry-specific selection
  const arrowSelection = selections[geometryKey] || [];
  const selectionSet = useMemo(() => new Set(arrowSelection), [arrowSelection]);

  // Individual queries for each attribute - enables perfect cross-component caching
  const { data: positionData, isFetching: isPositionFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, position],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [position as string], signal),
    enabled: !!roomId && !!clientId && typeof position === "string",
    placeholderData: keepPreviousData,
  });

  const { data: directionData, isFetching: isDirectionFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, direction],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [direction as string], signal),
    enabled: !!roomId && !!clientId && typeof direction === "string",
    placeholderData: keepPreviousData,
  });

  const { data: colorData, isFetching: isColorFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, color],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [color as string], signal),
    enabled: !!roomId && !!clientId && typeof color === "string" && shouldFetchAsFrameData(color as string),
    placeholderData: keepPreviousData,
  });

  const { data: radiusData, isFetching: isRadiusFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, radius],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [radius as string], signal),
    enabled: !!roomId && !!clientId && typeof radius === "string",
    placeholderData: keepPreviousData,
  });

  const { data: scaleData, isFetching: isScaleFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, scale],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [scale as string], signal),
    enabled: !!roomId && !!clientId && typeof scale === "string",
    placeholderData: keepPreviousData,
  });

  // Check if any enabled query is still fetching
  const isFetching =
    (typeof position === "string" && isPositionFetching) ||
    (typeof direction === "string" && isDirectionFetching) ||
    (typeof color === "string" && shouldFetchAsFrameData(color as string) && isColorFetching) ||
    (typeof radius === "string" && isRadiusFetching) ||
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

  useEffect(() => {
    if (isFetching) {
      return; // Wait for all enabled queries to complete
    }

    try {
      // --- Data Processing Step ---
      // Get fetched data or use static values
      const fetchedPosition = typeof position === 'string' ? positionData?.[position as string] : undefined;
      const finalCount = getInstanceCount(position, fetchedPosition);

      if (finalCount === 0) {
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      // Process all attributes using utility functions
      const finalPositions = processNumericAttribute(position, fetchedPosition, finalCount);

      const fetchedColor = typeof color === 'string' ? colorData?.[color as string] : undefined;
      const finalColors = processColorAttribute(color, fetchedColor, finalCount);

      const fetchedRadius = typeof radius === 'string' ? radiusData?.[radius as string] : undefined;
      const finalRadii = processNumericAttribute(radius, fetchedRadius, finalCount);

      const fetchedDirection = typeof direction === 'string' ? directionData?.[direction as string] : undefined;
      const finalDirections = processNumericAttribute(direction, fetchedDirection, finalCount);

      const fetchedScale = typeof scale === 'string' ? scaleData?.[scale as string] : undefined;
      const finalScales = processNumericAttribute(scale, fetchedScale, finalCount);

      // --- Validation Step ---
      const isDataValid = validateArrayLengths(
        {
          positions: finalPositions,
          directions: finalDirections,
          colors: finalColors,
          radii: finalRadii,
          scales: finalScales,
        },
        {
          positions: finalCount * 3,
          directions: finalCount * 3,
          colors: finalCount * 3,
          radii: finalCount,
          scales: finalCount,
        }
      );

      if (!isDataValid) {
        console.error("Arrow data is invalid or has inconsistent lengths.");
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      // --- Mesh Resizing Step ---
      if (instanceCount !== finalCount) {
        setInstanceCount(finalCount);
        return;
      }

      const mesh = instancedMeshRef.current;
      if (!mesh) return;

      // --- Mesh Instance Update Step ---
      for (let i = 0; i < finalCount; i++) {
        const i3 = i * 3;
        _vec3.set(finalPositions[i3], finalPositions[i3 + 1], finalPositions[i3 + 2]);
        _vec3_2.set(finalDirections[i3], finalDirections[i3 + 1], finalDirections[i3 + 2]);

        const arrowLength = _vec3_2.length() * finalScales[i];
        const arrowRadius = finalRadii[i];
        _vec3_3.set(arrowRadius, arrowLength, arrowRadius);

        // Avoid issues with zero-length vectors
        if (arrowLength > 1e-6) {
          _quat.setFromUnitVectors(_arrowUp, _vec3_2.normalize());
        } else {
          _quat.identity(); // No rotation
        }

        _matrix.compose(_vec3, _quat, _vec3_3);
        mesh.setMatrixAt(i, _matrix);

        if (selectionSet && selectionSet.has(i)) {
          _color.setRGB(SELECTION_COLOR[0], SELECTION_COLOR[1], SELECTION_COLOR[2]);
        } else {
          _color.setRGB(finalColors[i3], finalColors[i3 + 1], finalColors[i3 + 2]);
        }
        mesh.setColorAt(i, _color);
      }

      mesh.instanceMatrix.needsUpdate = true;
      if (mesh.instanceColor) {
        mesh.instanceColor.needsUpdate = true;
      }
    } catch (error) {
      console.error("Error processing Arrow data:", error);
      if (instanceCount !== 0) setInstanceCount(0);
    }
  }, [
    isFetching,
    positionData,
    directionData,
    colorData,
    radiusData,
    scaleData,
    position,
    direction,
    color,
    radius,
    scale,
    instanceCount,
    selectionSet,
  ]);

  // Convert instanced mesh to merged mesh for path tracing
  useEffect(() => {
    if (!pathtracingEnabled) {
      // Clean up merged mesh when pathtracing disabled
      if (mergedMeshRef.current) {
        disposeMesh(mergedMeshRef.current);
        mergedMeshRef.current = null;
      }
      return;
    }

    if (!instancedMeshRef.current || instanceCount === 0) return;

    // Dispose old merged mesh if it exists
    if (mergedMeshRef.current) {
      disposeMesh(mergedMeshRef.current);
    }

    // Convert instanced mesh to single merged mesh with vertex colors
    const mergedMesh = convertInstancedMeshToMerged(instancedMeshRef.current);
    mergedMeshRef.current = mergedMesh;

    // Request pathtracing update
    requestPathtracingUpdate();

    // Cleanup on unmount or when dependencies change
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
    // DO NOT depend on positionData/directionData/colorData/selections here!
    // That causes unnecessary rebuilds. instanceCount only changes AFTER mesh update completes.
  ]);

  // Create the base geometry only once
  const geometry = useMemo(createArrowMesh, []);

  const onClickHandler = useCallback((event: any) => {
    if (event.detail !== 1 || event.instanceId === undefined) return;
    event.stopPropagation();
    updateSelections(geometryKey, event.instanceId, event.shiftKey);
  }, [updateSelections, geometryKey]);

  if (!clientId || !roomId) {
    return null;
  }

  return (
    <group>
      {/* Instanced mesh - visible when NOT pathtracing */}
      {/* NOTE: Interactions (click) disabled when pathtracing enabled */}
      <instancedMesh
        key={instanceCount}
        ref={instancedMeshRef}
        args={[geometry, undefined, instanceCount]}
        visible={!pathtracingEnabled}
        onClick={!pathtracingEnabled ? onClickHandler : undefined}
      >
        {renderMaterial(material)}
      </instancedMesh>

      {/* Merged mesh - visible when pathtracing */}
      {pathtracingEnabled && mergedMeshRef.current && (
        <primitive object={mergedMeshRef.current} />
      )}
    </group>
  );
}
