import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { getFrames, createGeometry } from "../../myapi/client";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect, useCallback } from "react";
import { debounce } from "lodash";
import { useGeometryEditing } from "../../hooks/useGeometryEditing";
import { renderMaterial } from "./materials";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import {
  processPositionAttribute,
  processSize2D,
  processRotationAttribute,
  processColorData,
  getInstanceCount,
  validateArrayLengths,
  expandSharedColor,
  SELECTION_SCALE,
  HOVER_SCALE,
} from "../../utils/geometryData";
import { _vec3, _vec3_2, _vec3_3, _euler, _matrix, _matrix2, _quat, _quat2, _color } from "../../utils/threeObjectPools";
import { convertInstancedMeshToMerged, disposeMesh } from "../../utils/convertInstancedMesh";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";
import { useFrameKeys } from "../../hooks/useSchemas";

interface InteractionSettings {
  enabled: boolean;
  color: string;
  opacity: number;
}

interface PlaneData {
  position: string | number[][];
  size: string | number[][];
  color: string | string[]; // Dynamic ref or list of hex strings
  rotation: string | number[][];
  material: string;
  scale: number;
  opacity: number;
  double_sided: boolean;
  selecting: InteractionSettings;
  hovering: InteractionSettings;
  active?: boolean; // Whether geometry is active (can be disabled on critical errors)
}

export default function Plane({
  data,
  geometryKey,
  pathtracingEnabled = false
}: {
  data: PlaneData;
  geometryKey: string;
  pathtracingEnabled?: boolean;
}) {
  // Use individual selectors to prevent unnecessary re-renders
  const geometryDefaults = useAppStore((state) => state.geometryDefaults);

  // Merge with defaults from Pydantic (single source of truth)
  const fullData = getGeometryWithDefaults<PlaneData>(data, "Plane", geometryDefaults);

  const {
    position: positionProp,
    size: sizeProp,
    color: colorProp,
    rotation: rotationProp,
    material,
    scale,
    double_sided,
    selecting,
    hovering,
    opacity,
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
    if (typeof sizeProp === "string") keys.push(sizeProp);
    return keys;
  }, [positionProp, sizeProp]);

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
  const planeSelection = selections[geometryKey] || [];
  const selectionSet = useMemo(() => new Set(planeSelection), [planeSelection]);
  const selectedIndices = useMemo(() => Array.from(selectionSet), [selectionSet]);
  const validSelectedIndices = useMemo(
    () => selectedIndices.filter((id) => id < instanceCount),
    [selectedIndices, instanceCount]
  );

  const planeScale = scale || 1.0;

  // Individual queries for each attribute
  const { data: positionData, isFetching: isPositionFetching, isError: isPositionError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, positionProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [positionProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof positionProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: sizeData, isFetching: isSizeFetching, isError: isSizeError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, sizeProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [sizeProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof sizeProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: colorData, isFetching: isColorFetching, isError: isColorError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, colorProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [colorProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string),
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: rotationData, isFetching: isRotationFetching, isError: isRotationError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, rotationProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [rotationProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof rotationProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  // Check if any enabled query is still fetching
  const isFetching =
    (typeof positionProp === "string" && isPositionFetching) ||
    (typeof sizeProp === "string" && isSizeFetching) ||
    (typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string) && isColorFetching) ||
    (typeof rotationProp === "string" && isRotationFetching);

  // Check if any query has errored - treat as data unavailable
  const hasQueryError = useMemo(
    () =>
      (typeof positionProp === "string" && isPositionError) ||
      (typeof sizeProp === "string" && isSizeError) ||
      (typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string) && isColorError) ||
      (typeof rotationProp === "string" && isRotationError),
    [positionProp, isPositionError, sizeProp, isSizeError, colorProp, isColorError, rotationProp, isRotationError]
  );

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
    "Plane",
    fullData
  );

  // Persist position changes to server (debounced)
  // Watch the geometry's position in Zustand store and persist when it changes locally
  const persistPositions = useCallback(async () => {
    if (!roomId) return;

    // Get current geometry from Zustand store
    const currentGeometry = geometries[geometryKey];
    if (!currentGeometry || !currentGeometry.data) return;

    const currentPosition = currentGeometry.data.position;

    // Only persist if position is static (number[][])
    if (!Array.isArray(currentPosition) || currentPosition.length === 0 || !Array.isArray(currentPosition[0])) {
      return;
    }

    try {
      await createGeometry(roomId, geometryKey, "Plane", currentGeometry.data, lock?.token);
    } catch (error: any) {
      console.error(`[Plane] Failed to persist ${geometryKey}:`, error);
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

  // Consolidated data processing and mesh update
  useEffect(() => {
    if (isFetching) {
      return; // Wait for all enabled queries to complete
    }

    // If queries have errored, continue with fallback to static data (this is normal when data doesn't exist)
    // No logging needed as this is expected behavior

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
      const colorHexArray = processColorData(colorProp, fetchedColor, finalCount);

      const fetchedSize = typeof sizeProp === 'string' ? sizeData?.[sizeProp as string] : undefined;
      const finalSizes = processSize2D(sizeProp, fetchedSize, finalCount);

      const fetchedRotation = typeof rotationProp === 'string' ? rotationData?.[rotationProp as string] : undefined;
      const finalRotations = processRotationAttribute(rotationProp, fetchedRotation, finalCount);

      // Handle shared color (single color for all instances)
      const finalColorHex = expandSharedColor(colorHexArray, finalCount);

      // --- Validation Step ---
      const isDataValid = validateArrayLengths(
        { positions: finalPositions, sizes: finalSizes, rotations: finalRotations },
        { positions: finalCount * 3, sizes: finalCount * 2, rotations: finalCount * 3 }
      ) && (finalColorHex.length === finalCount);

      if (!isDataValid) {
        console.error("Plane data is invalid or has inconsistent lengths.");
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
        const i2 = i * 2;
        _vec3.set(finalPositions[i3], finalPositions[i3 + 1], finalPositions[i3 + 2]);
        _euler.set(finalRotations[i3], finalRotations[i3 + 1], finalRotations[i3 + 2]);
        const width = finalSizes[i2] * planeScale;
        const height = finalSizes[i2 + 1] * planeScale;
        _quat.setFromEuler(_euler);
        _vec3_2.set(width, height, 1);
        _matrix.compose(_vec3, _quat, _vec3_2);
        mainMesh.setMatrixAt(i, _matrix);

        // Set color directly from hex string (THREE.Color.set() accepts hex)
        _color.set(finalColorHex[i]);
        mainMesh.setColorAt(i, _color);
      }

      mainMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
      mainMesh.instanceMatrix.needsUpdate = true;
      if (mainMesh.instanceColor) mainMesh.instanceColor.needsUpdate = true;

      // Update bounding box to prevent frustum culling issues
      mainMesh.computeBoundingBox();
      mainMesh.computeBoundingSphere();

      // --- Selection Mesh Update (uses thin box geometry for visible outline) ---
      if (selecting.enabled && selectionMeshRef.current) {
        const selectionMesh = selectionMeshRef.current;
        validSelectedIndices.forEach((id, index) => {
          if (id >= finalCount) return;
          const i3 = id * 3;
          const i2 = id * 2;
          _vec3.set(finalPositions[i3], finalPositions[i3 + 1], finalPositions[i3 + 2]);
          _euler.set(finalRotations[i3], finalRotations[i3 + 1], finalRotations[i3 + 2]);
          const width = finalSizes[i2] * planeScale * SELECTION_SCALE;
          const height = finalSizes[i2 + 1] * planeScale * SELECTION_SCALE;
          // Use thin depth (0.01) for visible outline effect
          _quat.setFromEuler(_euler);
          _vec3_2.set(width, height, 0.01);
          _matrix.compose(_vec3, _quat, _vec3_2);
          selectionMesh.setMatrixAt(index, _matrix);
        });
        selectionMesh.instanceMatrix.needsUpdate = true;

        // Update bounding box for selection mesh
        selectionMesh.computeBoundingBox();
        selectionMesh.computeBoundingSphere();
      }

    } catch (error) {
      console.error("Error processing Plane data:", error);
      if (instanceCount !== 0) setInstanceCount(0);
    }
  }, [
    data, // Add data to dependencies to ensure updates trigger
    isFetching,
    hasQueryError,
    positionData,
    sizeData,
    colorData,
    rotationData,
    positionProp,
    sizeProp,
    colorProp,
    rotationProp,
    instanceCount,
    planeScale,
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

      // Apply hover scale (with thin depth for outline)
      hoverMesh.position.copy(_vec3);
      hoverMesh.quaternion.copy(_quat2);
      hoverMesh.scale.set(
        _vec3_2.x * HOVER_SCALE,
        _vec3_2.y * HOVER_SCALE,
        0.01 // Thin depth for outline
      );
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

  // Shared geometries
  const mainGeometry = useMemo(() => new THREE.PlaneGeometry(1, 1), []);
  // Thin box geometry for selection/hover (since PlaneGeometry can't be visibly scaled for outline)
  const outlineGeometry = useMemo(() => new THREE.BoxGeometry(1, 1, 0.01), []);

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

  if (!userName || !roomId) return null;

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
        args={[undefined, undefined, instanceCount]}
        visible={!pathtracingEnabled}
        onClick={!pathtracingEnabled && selecting.enabled ? onClickHandler : undefined}
        onPointerEnter={!pathtracingEnabled && hovering?.enabled ? onPointerEnterHandler : undefined}
        onPointerMove={!pathtracingEnabled && hovering?.enabled ? onPointerMoveHandler : undefined}
        onPointerOut={!pathtracingEnabled && hovering?.enabled ? onPointerOutHandler : undefined}
      >
        <primitive object={mainGeometry} attach="geometry" />
        {renderMaterial(material, opacity, undefined, double_sided ? THREE.DoubleSide : THREE.FrontSide)}
      </instancedMesh>

      {/* Selection mesh - thin box geometry for visible outline */}
      {!pathtracingEnabled && selecting.enabled && (
        <instancedMesh
          key={`selection-${validSelectedIndices.length}`}
          ref={selectionMeshRef}
          args={[undefined, undefined, validSelectedIndices.length]}
        >
          <primitive object={outlineGeometry} attach="geometry" />
          <meshBasicMaterial
            side={double_sided ? THREE.DoubleSide : THREE.FrontSide}
            transparent
            opacity={selecting.opacity}
            color={selecting.color}
          />
        </instancedMesh>
      )}

      {/* Hover mesh - thin box geometry for visible outline */}
      {!pathtracingEnabled && hovering?.enabled && (
        <mesh ref={hoverMeshRef} visible={false}>
          <primitive object={outlineGeometry} attach="geometry" />
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
