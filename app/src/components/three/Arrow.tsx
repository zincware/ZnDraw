import * as THREE from "three";
import { useAppStore } from "../../store";
import { useRef, useMemo, useEffect, useState, useCallback } from "react";
import { useGeometryEditing } from "../../hooks/useGeometryEditing";
import { BufferGeometryUtils } from "three/examples/jsm/Addons.js";
import { renderMaterial } from "./materials";
import { getFrames, createGeometry } from "../../myapi/client";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { debounce } from "lodash";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import {
  type SizeProp,
  processNumericAttribute,
  processColorData,
  getInstanceCount,
  validateArrayLengths,
  expandSharedColor,
  SELECTION_SCALE,
  HOVER_SCALE,
} from "../../utils/geometryData";
import { _vec3, _vec3_2, _vec3_3, _quat, _quat2, _matrix, _matrix2, _color } from "../../utils/threeObjectPools";
import { convertInstancedMeshToMerged, disposeMesh } from "../../utils/convertInstancedMesh";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";
import { useFrameKeys } from "../../hooks/useSchemas";

interface InteractionSettings {
  enabled: boolean;
  color: string;
  opacity: number;
}

interface ArrowData {
  position: string | number[][];
  direction: string | number[][];
  color: string | string[]; // Dynamic ref or list of hex strings
  radius: SizeProp;
  scale: SizeProp;
  material: string;
  opacity: number;
  selecting: InteractionSettings;
  hovering: InteractionSettings;
  active?: boolean; // Whether geometry is active (can be disabled on critical errors)
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
  data: ArrowData;
  geometryKey: string;
  pathtracingEnabled?: boolean;
}) {
  const geometryDefaults = useAppStore((state) => state.geometryDefaults);

  // Merge with defaults from Pydantic (single source of truth)
  const fullData = getGeometryWithDefaults<ArrowData>(data, "Arrow", geometryDefaults);
  const {
    position: positionProp,
    direction: directionProp,
    color,
    radius,
    scale,
    material,
    opacity,
    selecting,
    hovering
  } = fullData;

  const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const hoverMeshRef = useRef<THREE.Mesh | null>(null);
  const mergedMeshRef = useRef<THREE.Mesh | null>(null);
  const [instanceCount, setInstanceCount] = useState(0);

  // Use individual selectors to prevent unnecessary re-renders
  const currentFrame = useAppStore((state) => state.currentFrame);
  const frameCount = useAppStore((state) => state.frameCount);
  const roomId = useAppStore((state) => state.roomId);
  const userName = useAppStore((state) => state.userName);
  const lock = useAppStore((state) => state.lock);
  const selections = useAppStore((state) => state.selections);
  const updateSelections = useAppStore((state) => state.updateSelections);
  const hoveredGeometryInstance = useAppStore((state) => state.hoveredGeometryInstance);
  const setHoveredGeometryInstance = useAppStore((state) => state.setHoveredGeometryInstance);
  const setDrawingPointerPosition = useAppStore((state) => state.setDrawingPointerPosition);
  const mode = useAppStore((state) => state.mode);
  const setDrawingIsValid = useAppStore((state) => state.setDrawingIsValid);
  const setGeometryFetching = useAppStore((state) => state.setGeometryFetching);
  const removeGeometryFetching = useAppStore((state) => state.removeGeometryFetching);
  const requestPathtracingUpdate = useAppStore((state) => state.requestPathtracingUpdate);
  const showSnackbar = useAppStore((state) => state.showSnackbar);
  const geometries = useAppStore((state) => state.geometries);
  const geometryUpdateSources = useAppStore((state) => state.geometryUpdateSources);

  // Fetch frame keys to check if required data is available
  const { data: frameKeysData, isLoading: isLoadingKeys } = useFrameKeys(roomId!, currentFrame);

  // Check if required keys are available for this geometry
  const requiredKeys = useMemo(() => {
    const keys: string[] = [];
    if (typeof positionProp === "string") keys.push(positionProp);
    if (typeof directionProp === "string") keys.push(directionProp);
    return keys;
  }, [positionProp, directionProp]);

  const hasRequiredKeys = useMemo(() => {
    // While loading keys, assume we have them (keep previous frame rendered)
    if (isLoadingKeys) return true;
    // If no keys data yet, assume we have them (keep previous frame)
    if (!frameKeysData?.keys) return true;
    // Only when we have keys data, check if required keys are available
    const availableKeys = new Set(frameKeysData.keys);
    return requiredKeys.every(key => availableKeys.has(key));
  }, [frameKeysData, requiredKeys, isLoadingKeys]);

  // Use geometry-specific selection
  const arrowSelection = selections[geometryKey] || [];
  const selectionSet = useMemo(() => new Set(arrowSelection), [arrowSelection]);
  const selectedIndices = useMemo(() => Array.from(selectionSet), [selectionSet]);
  const validSelectedIndices = useMemo(
    () => selectedIndices.filter((id) => id < instanceCount),
    [selectedIndices, instanceCount]
  );

  // Individual queries for each attribute - enables perfect cross-component caching
  const { data: positionData, isFetching: isPositionFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, positionProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [positionProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof positionProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: directionData, isFetching: isDirectionFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, directionProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [directionProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof directionProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: colorData, isFetching: isColorFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, color],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [color as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof color === "string" && shouldFetchAsFrameData(color as string),
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: radiusData, isFetching: isRadiusFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, radius],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [radius as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof radius === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: scaleData, isFetching: isScaleFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, scale],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [scale as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof scale === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  // Check if any enabled query is still fetching
  const isFetching =
    (typeof positionProp === "string" && isPositionFetching) ||
    (typeof directionProp === "string" && isDirectionFetching) ||
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

  // Handle geometry editing with transform controls
  const finalPositionData = typeof positionProp === "string" ? positionData?.[positionProp] : positionProp;
  useGeometryEditing(
    geometryKey,
    finalPositionData,
    selectedIndices,
    "Arrow",
    fullData
  );

  // Persistence callback - persists position changes to server
  const persistPositions = useCallback(async () => {
    if (!roomId) return;

    const currentGeometry = geometries[geometryKey];
    if (!currentGeometry || !currentGeometry.data) return;

    const currentPosition = currentGeometry.data.position;

    // Only persist if position is static (number[][])
    if (!Array.isArray(currentPosition) || currentPosition.length === 0 || !Array.isArray(currentPosition[0])) {
      return;
    }

    try {
      await createGeometry(roomId, geometryKey, "Arrow", currentGeometry.data, lock?.token);
    } catch (error: any) {
      console.error(`[Arrow] Failed to persist ${geometryKey}:`, error);
      // Snackbar is shown automatically by withAutoLock for lock failures
    }
  }, [roomId, geometryKey, geometries, lock, showSnackbar]);

  // Memoize debounced persist function to avoid recreation on every render
  const debouncedPersist = useMemo(
    () => debounce(persistPositions, 500),
    [persistPositions]
  );

  // Cleanup debounce on unmount
  useEffect(() => {
    return () => {
      debouncedPersist.cancel();
    };
  }, [debouncedPersist]);

  // Watch position changes and persist - only if source is 'local'
  useEffect(() => {
    const currentGeometry = geometries[geometryKey];
    if (!currentGeometry) return;

    const currentPosition = currentGeometry.data?.position;
    if (!currentPosition) return;

    // Only persist if position is static
    if (!Array.isArray(currentPosition) || currentPosition.length === 0 || !Array.isArray(currentPosition[0])) {
      return;
    }

    // Only persist if update source is 'local' (not from server)
    const updateSource = geometryUpdateSources[geometryKey];
    if (updateSource !== 'local') {
      return;
    }

    debouncedPersist();
  }, [geometries[geometryKey]?.data?.position, geometryUpdateSources[geometryKey], debouncedPersist, geometryKey]);


  useEffect(() => {
    if (isFetching) {
      return; // Wait for all enabled queries to complete
    }

    // If queries have errored, continue with fallback to static data (this is normal when data doesn't exist)
    // No logging needed as this is expected behavior

    try {
      // --- Data Processing Step ---
      // Get fetched data or use static values
      const fetchedPosition = typeof positionProp === 'string' ? positionData?.[positionProp as string] : undefined;
      const finalCount = getInstanceCount(positionProp, fetchedPosition);

      if (finalCount === 0) {
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      // Process all attributes using utility functions
      const finalPositions = processNumericAttribute(positionProp, fetchedPosition, finalCount);

      const fetchedColor = typeof color === 'string' ? colorData?.[color as string] : undefined;
      const colorHexArray = processColorData(color, fetchedColor, finalCount);

      const fetchedRadius = typeof radius === 'string' ? radiusData?.[radius as string] : undefined;
      const finalRadii = processNumericAttribute(radius, fetchedRadius, finalCount);

      const fetchedDirection = typeof directionProp === 'string' ? directionData?.[directionProp as string] : undefined;
      const finalDirections = processNumericAttribute(directionProp, fetchedDirection, finalCount);

      const fetchedScale = typeof scale === 'string' ? scaleData?.[scale as string] : undefined;
      const finalScales = processNumericAttribute(scale, fetchedScale, finalCount);

      // Handle shared color (single color for all instances)
      const finalColorHex = expandSharedColor(colorHexArray, finalCount);

      // --- Validation Step ---
      const isDataValid = validateArrayLengths(
        {
          positions: finalPositions,
          directions: finalDirections,
          radii: finalRadii,
          scales: finalScales,
        },
        {
          positions: finalCount * 3,
          directions: finalCount * 3,
          radii: finalCount,
          scales: finalCount,
        }
      ) && (finalColorHex.length === finalCount);

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

      // --- Main Mesh Instance Update ---
      const mainMesh = mainMeshRef.current;
      if (!mainMesh) return;

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
        mainMesh.setMatrixAt(i, _matrix);

        // Set color directly from hex string (THREE.Color.set() accepts hex)
        _color.set(finalColorHex[i]);
        mainMesh.setColorAt(i, _color);
      }

      mainMesh.instanceMatrix.needsUpdate = true;
      if (mainMesh.instanceColor) {
        mainMesh.instanceColor.needsUpdate = true;
      }

      // Update bounding box to prevent frustum culling issues
      mainMesh.computeBoundingBox();
      mainMesh.computeBoundingSphere();

      // --- Selection Mesh Update ---
      if (selecting.enabled && selectionMeshRef.current) {
        const selectionMesh = selectionMeshRef.current;
        validSelectedIndices.forEach((id, index) => {
          if (id >= finalCount) return;
          const i3 = id * 3;
          _vec3.set(finalPositions[i3], finalPositions[i3 + 1], finalPositions[i3 + 2]);
          _vec3_2.set(finalDirections[i3], finalDirections[i3 + 1], finalDirections[i3 + 2]);

          const arrowLength = _vec3_2.length() * finalScales[id];
          const arrowRadius = finalRadii[id] * SELECTION_SCALE;
          _vec3_3.set(arrowRadius, arrowLength * SELECTION_SCALE, arrowRadius);

          // Avoid issues with zero-length vectors
          if (arrowLength > 1e-6) {
            _quat.setFromUnitVectors(_arrowUp, _vec3_2.normalize());
          } else {
            _quat.identity();
          }

          _matrix.compose(_vec3, _quat, _vec3_3);
          selectionMesh.setMatrixAt(index, _matrix);
        });
        selectionMesh.instanceMatrix.needsUpdate = true;

        // Update bounding box for selection mesh
        selectionMesh.computeBoundingBox();
        selectionMesh.computeBoundingSphere();
      }
    } catch (error) {
      console.error("Error processing Arrow data:", error);
      if (instanceCount !== 0) setInstanceCount(0);
    }
  }, [
    data, // Add data to dependencies to ensure updates trigger
    isFetching,
    positionData,
    directionData,
    colorData,
    radiusData,
    scaleData,
    positionProp,
    directionProp,
    color,
    radius,
    scale,
    instanceCount,
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

      // Get transform from main mesh using pooled objects
      mainMesh.getMatrixAt(hoveredGeometryInstance.instanceId, _matrix2);
      _matrix2.decompose(_vec3, _quat2, _vec3_2);

      // Apply hover scale
      hoverMesh.position.copy(_vec3);
      hoverMesh.quaternion.copy(_quat2);
      hoverMesh.scale.set(_vec3_2.x * HOVER_SCALE, _vec3_2.y, _vec3_2.z * HOVER_SCALE);
    } else {
      hoverMesh.visible = false;
    }
  }, [hoveredGeometryInstance, instanceCount, hovering, geometryKey]);

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

    if (!mainMeshRef.current || instanceCount === 0) return;

    // Dispose old merged mesh if it exists
    if (mergedMeshRef.current) {
      disposeMesh(mergedMeshRef.current);
    }

    // Convert instanced mesh to single merged mesh with vertex colors
    const mergedMesh = convertInstancedMeshToMerged(mainMeshRef.current);
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
  ]);

  // Create the base geometry only once
  const geometry = useMemo(createArrowMesh, []);

  const onClickHandler = useCallback((event: any) => {
    if (event.detail !== 1 || event.instanceId === undefined) return;
    event.stopPropagation();
    updateSelections(geometryKey, event.instanceId, event.shiftKey);
  }, [updateSelections, geometryKey]);

  const onPointerMoveHandler = useCallback((event: any) => {
    if (event.instanceId === undefined) return;
    event.stopPropagation();
    if (mode === 'drawing') {
      setDrawingPointerPosition(event.point);
    }
  }, [mode, setDrawingPointerPosition]);

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

  if (!userName || !roomId) {
    return null;
  }

  // Don't render if geometry is disabled OR if required keys are not available
  if (fullData.active === false || !hasRequiredKeys) {
    return null;
  }

  return (
    <group>
      {/* Main instanced mesh - visible when NOT pathtracing */}
      <instancedMesh
        key={instanceCount}
        ref={mainMeshRef}
        args={[geometry, undefined, instanceCount]}
        visible={!pathtracingEnabled}
        onClick={!pathtracingEnabled && selecting.enabled ? onClickHandler : undefined}
        onPointerEnter={!pathtracingEnabled && hovering?.enabled ? onPointerEnterHandler : undefined}
        onPointerMove={!pathtracingEnabled && hovering?.enabled ? onPointerMoveHandler : undefined}
        onPointerOut={!pathtracingEnabled && hovering?.enabled ? onPointerOutHandler : undefined}
      >
        {renderMaterial(material, opacity)}
      </instancedMesh>

      {/* Selection mesh - only when NOT pathtracing */}
      {!pathtracingEnabled && selecting.enabled && (
        <instancedMesh
          key={`selection-${validSelectedIndices.length}`}
          ref={selectionMeshRef}
          args={[geometry, undefined, validSelectedIndices.length]}
        >
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
          <primitive object={geometry} attach="geometry" />
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
