import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { getFrames, updateGeometryActive, createGeometry } from "../../myapi/client";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect, useCallback } from "react";
import { useGeometryEditing } from "../../hooks/useGeometryEditing";
import { renderMaterial } from "./materials";
import { debounce } from "lodash";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import {
  processNumericAttribute,
  processColorData,
  getInstanceCount,
  validateArrayLengths,
  expandSharedColor,
  SELECTION_SCALE,
  HOVER_SCALE,
} from "../../utils/geometryData";
import { _vec3, _vec3_2, _matrix, _matrix2, _quat2, _color } from "../../utils/threeObjectPools";
import { convertInstancedMeshToMerged, disposeMesh } from "../../utils/convertInstancedMesh";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";
import { isTransform, evaluateTransform, getTransformSources, type Transform } from "../../utils/transformProcessor";

interface InteractionSettings {
  enabled: boolean;
  color: string;
  opacity: number;
}

interface SphereData {
  position: string | number[][] | Transform;
  color: string | string[] | Transform; // Dynamic ref or list of hex strings or transform
  radius: string | number[] | number | Transform;
  material: string;
  resolution: number;
  scale: number;
  opacity: number;
  selecting: InteractionSettings;
  hovering: InteractionSettings;
  active?: boolean; // Whether geometry is active (can be disabled on critical errors)
}

// Reusable THREE objects imported from threeObjectPools

export default function Sphere({
  data,
  geometryKey,
  pathtracingEnabled = false
}: {
  data: SphereData;
  geometryKey: string;
  pathtracingEnabled?: boolean;
}) {
  const { geometryDefaults, currentFrame, frameCount, roomId, userName, selections, updateSelections, setDrawingPointerPosition, isDrawing, setDrawingIsValid, setGeometryFetching, removeGeometryFetching, hoveredGeometryInstance, setHoveredGeometryInstance, setParticleCount, requestPathtracingUpdate, updateGeometry, showSnackbar, geometries, geometryUpdateSources } = useAppStore();

  // Merge with defaults from Pydantic (single source of truth)
  const fullData = getGeometryWithDefaults<SphereData>(data, "Sphere", geometryDefaults);

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
  } = fullData;

  const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const hoverMeshRef = useRef<THREE.Mesh | null>(null);
  const mergedMeshRef = useRef<THREE.Mesh | null>(null);
  const [instanceCount, setInstanceCount] = useState(0);
  const hasDisabledGeometryRef = useRef(false);


  // Use geometry-specific selection
  const particleSelection = selections[geometryKey] || [];
  const selectionSet = useMemo(() => new Set(particleSelection), [particleSelection]);
  const selectedIndices = useMemo(() => Array.from(selectionSet), [selectionSet]);
  const validSelectedIndices = useMemo(
    () => selectedIndices.filter((id) => id < instanceCount),
    [selectedIndices, instanceCount]
  );

  const particleResolution = resolution || 8;
  const particleScale = scale || 1.0;

  // Check if properties are transforms and get required sources
  const positionIsTransform = isTransform(positionProp);
  const colorIsTransform = isTransform(colorProp);
  const radiusIsTransform = isTransform(radiusProp);

  const positionTransformSources = positionIsTransform ? getTransformSources(positionProp as Transform) : [];
  const colorTransformSources = colorIsTransform ? getTransformSources(colorProp as Transform) : [];
  const radiusTransformSources = radiusIsTransform ? getTransformSources(radiusProp as Transform) : [];

  // Collect all unique keys needed for transforms
  const allTransformKeys = useMemo(() => {
    const keys = new Set<string>();
    [...positionTransformSources, ...colorTransformSources, ...radiusTransformSources].forEach(k => keys.add(k));
    return Array.from(keys);
  }, [positionTransformSources, colorTransformSources, radiusTransformSources]);

  // Fetch all transform source data in a single query
  const { data: transformData, isFetching: isTransformFetching, isError: isTransformError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, ...allTransformKeys],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, allTransformKeys, signal),
    enabled: !!roomId && !!userName && frameCount > 0 && allTransformKeys.length > 0,
    placeholderData: keepPreviousData,
    retry: false,
  });

  // Individual queries for each attribute - enables perfect cross-component caching
  const { data: positionData, isFetching: isPositionFetching, isError: isPositionError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, positionProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [positionProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof positionProp === "string" && !positionIsTransform,
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: colorData, isFetching: isColorFetching, isError: isColorError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, colorProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [colorProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof colorProp === "string" && !colorIsTransform && shouldFetchAsFrameData(colorProp as string),
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: radiusData, isFetching: isRadiusFetching, isError: isRadiusError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, radiusProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [radiusProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof radiusProp === "string" && !radiusIsTransform,
    placeholderData: keepPreviousData,
    retry: false,
  });

  // Check if any enabled query is still fetching
  const isFetching =
    (allTransformKeys.length > 0 && isTransformFetching) ||
    (typeof positionProp === "string" && !positionIsTransform && isPositionFetching) ||
    (typeof colorProp === "string" && !colorIsTransform && shouldFetchAsFrameData(colorProp as string) && isColorFetching) ||
    (typeof radiusProp === "string" && !radiusIsTransform && isRadiusFetching);

  // Check if any query has errored - treat as data unavailable
  const hasQueryError = useMemo(
    () =>
      (allTransformKeys.length > 0 && isTransformError) ||
      (typeof positionProp === "string" && !positionIsTransform && isPositionError) ||
      (typeof colorProp === "string" && !colorIsTransform && shouldFetchAsFrameData(colorProp as string) && isColorError) ||
      (typeof radiusProp === "string" && !radiusIsTransform && isRadiusError),
    [allTransformKeys.length, isTransformError, positionProp, positionIsTransform, isPositionError, colorProp, colorIsTransform, isColorError, radiusProp, radiusIsTransform, isRadiusError]
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
    "Sphere",
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
      await createGeometry(roomId, geometryKey, "Sphere", currentGeometry.data);
    } catch (error) {
      console.error(`[Particles] Failed to persist ${geometryKey}:`, error);
      showSnackbar(`Failed to save positions for ${geometryKey}`, "error");
    }
  }, [roomId, geometryKey, geometries, showSnackbar]);

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

  // Detect critical fetch failures and disable geometry
  useEffect(() => {
    if (!roomId || !userName || hasDisabledGeometryRef.current || isFetching || frameCount === 0) {
      return;
    }

    // Critical error for Particles/Sphere: position query failed
    const hasCriticalError =
      (positionIsTransform && isTransformError) ||
      (typeof positionProp === "string" && !positionIsTransform && isPositionError);

    if (hasCriticalError) {
      if (process.env.NODE_ENV !== 'production') {
        console.warn(
          `Particles geometry "${geometryKey}": Critical data fetch failed. Disabling geometry.`,
          { positionError: isPositionError, transformError: isTransformError }
        );
      }
      hasDisabledGeometryRef.current = true;

      // Optimistically update local state immediately
      const updatedGeometry = {
        type: "Sphere",
        data: { ...data, active: false }
      };
      updateGeometry(geometryKey, updatedGeometry);

      // Show snackbar notification
      showSnackbar(`Geometry "${geometryKey}" disabled - data fetch failed`, "warning");

      // Then update server (server will skip emitting back to this client)
      updateGeometryActive(roomId, geometryKey, "Sphere", false)
        .then(() => {
          if (process.env.NODE_ENV !== 'production') {
            console.info(`Particles geometry "${geometryKey}" disabled successfully on server.`);
          }
        })
        .catch((error) => {
          console.error(`Failed to disable Particles geometry "${geometryKey}" on server:`, error);
          // Rollback optimistic update on error
          const rollbackGeometry = {
            type: "Sphere",
            data: { ...data }
          };
          updateGeometry(geometryKey, rollbackGeometry);
          hasDisabledGeometryRef.current = false;
        });
    }
  }, [
    roomId,
    userName,
    geometryKey,
    frameCount,
    isFetching,
    positionProp,
    positionIsTransform,
    isPositionError,
    isTransformError,
    updateGeometry,
    showSnackbar,
  ]);

  // Consolidated data processing and mesh update
  useEffect(() => {
    if (isFetching) {
      return; // Wait for all enabled queries to complete
    }

    // If queries have errored, continue with fallback to static data (this is normal when data doesn't exist)
    // No logging needed as this is expected behavior

    try {
      // --- Data Processing Step ---

      // If using transforms, wait for transform data to be loaded
      if ((positionIsTransform || radiusIsTransform || colorIsTransform) && !transformData) {
        // Transform data not yet loaded - skip rendering
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      // Evaluate transforms if needed
      let fetchedPosition: any;
      if (positionIsTransform && transformData) {
        fetchedPosition = evaluateTransform(positionProp as Transform, transformData);
      } else if (typeof positionProp === 'string') {
        fetchedPosition = positionData?.[positionProp as string];
      } else {
        fetchedPosition = undefined;
      }

      // Evaluate radius transform if needed
      let fetchedRadius: any;
      if (radiusIsTransform && transformData) {
        fetchedRadius = evaluateTransform(radiusProp as Transform, transformData);
      } else if (typeof radiusProp === 'string') {
        fetchedRadius = radiusData?.[radiusProp as string];
      } else {
        fetchedRadius = undefined;
      }

      // Evaluate color transform if needed
      let fetchedColor: any;
      if (colorIsTransform && transformData) {
        const colorFloats = evaluateTransform(colorProp as Transform, transformData);
        // Transform returns Float32Array of filtered values, but colors in zndraw are hex strings
        // For now, transforms on colors are not supported - fall back to static
        fetchedColor = undefined;
      } else if (typeof colorProp === 'string') {
        fetchedColor = colorData?.[colorProp as string];
      } else {
        fetchedColor = undefined;
      }

      // Calculate instance counts from each attribute
      const positionCount = getInstanceCount(positionProp, fetchedPosition);
      const radiusCount = fetchedRadius ? fetchedRadius.length : 0;

      // If ANY transform returns empty data, the geometry should render 0 instances
      // This handles cases where one transform succeeds but another fails/returns empty
      if ((positionIsTransform && positionCount === 0) ||
          (radiusIsTransform && radiusCount === 0)) {
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      const finalCount = positionCount;

      if (finalCount === 0) {
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      // Process all attributes
      const finalPositions = processNumericAttribute(positionProp, fetchedPosition, finalCount);
      const colorHexArray = processColorData(colorProp, fetchedColor, finalCount);
      const finalRadii = processNumericAttribute(radiusProp, fetchedRadius, finalCount);

      // Handle shared color (single color for all instances)
      const finalColorHex = expandSharedColor(colorHexArray, finalCount);

      // --- Validation Step ---
      const isDataValid = validateArrayLengths(
        { positions: finalPositions, radii: finalRadii },
        { positions: finalCount * 3, radii: finalCount }
      ) && (finalColorHex.length === finalCount);

      if (!isDataValid) {
        console.error("Sphere/Particles data is invalid or has inconsistent lengths.");
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      // --- Mesh Resizing Step ---
      if (instanceCount !== finalCount) {
        setInstanceCount(finalCount);
        setParticleCount(finalCount);
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

        // Update bounding box for selection mesh
        selectionMesh.computeBoundingBox();
        selectionMesh.computeBoundingSphere();
      }

    } catch (error) {
      console.error("Error processing Sphere/Particles data:", error);
      if (instanceCount !== 0) setInstanceCount(0);
    }
  }, [
    data, // Add data to dependencies to ensure updates trigger
    isFetching,
    hasQueryError,
    positionData,
    colorData,
    radiusData,
    transformData,
    positionProp,
    colorProp,
    radiusProp,
    positionIsTransform,
    radiusIsTransform,
    colorIsTransform,
    instanceCount,
    particleScale,
    validSelectedIndices,
    selecting,
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
      hoverMesh.scale.copy(_vec3_2).multiplyScalar(HOVER_SCALE);
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
    // DO NOT depend on positionData/colorData/radiusData/selections/hover here!
    // That causes unnecessary rebuilds. instanceCount only changes AFTER mesh update completes.
  ]);

  // Shared geometry for all particles (both instanced and merged meshes)
  const mainGeometry = useMemo(() => {
    return new THREE.SphereGeometry(1, particleResolution, particleResolution);
  }, [particleResolution]);

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
  }, [setDrawingIsValid, setHoveredGeometryInstance, geometryKey]);

  const onPointerOutHandler = useCallback(() => {
    setHoveredGeometryInstance(null, null);
    setDrawingIsValid(false);
  }, [setDrawingIsValid, setHoveredGeometryInstance]);

  if (!userName || !roomId) return null;

  // Don't render if geometry is disabled
  if (fullData.active === false) {
    return null;
  }

  return (
    <group>
      {/* Main instanced mesh - visible when NOT pathtracing */}
      {/* NOTE: Interactions (click, hover) disabled when pathtracing enabled */}
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
        {renderMaterial(material, opacity)}
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