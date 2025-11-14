import { useMemo, useRef, useState, useEffect, useCallback } from "react";
import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { useTheme } from "@mui/material/styles";
import { Line } from "@react-three/drei";
import { useAppStore } from "../../store";
import { getFrames, createGeometry } from "../../myapi/client";
import { debounce } from "lodash";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";
import { useFrameKeys } from "../../hooks/useSchemas";
import { useGeometryEditing } from "../../hooks/useGeometryEditing";
import {
  processPositionAttribute,
  getInstanceCount,
  validateArrayLengths,
  SELECTION_SCALE,
  HOVER_SCALE,
} from "../../utils/geometryData";
import { _vec3, _vec3_2, _matrix, _matrix2, _quat2 } from "../../utils/threeObjectPools";

interface InteractionSettings {
  enabled: boolean;
  color: string;
  opacity: number;
}

interface MarkerData {
  enabled: boolean;
  size: number;
  color: string;
  opacity: number;
  selecting?: InteractionSettings;
  hovering?: InteractionSettings;
}

interface CurveData {
  position: string | number[][];
  color: string;
  material: "LineBasicMaterial" | "LineDashedMaterial";
  divisions: number;
  variant: "CatmullRomCurve3";
  thickness: number;
  marker: MarkerData;
  virtual_marker: MarkerData;
  active?: boolean; // Whether geometry is active (can be disabled on critical errors)
}

/**
 * A component to render a 3D curve using @react-three/drei,
 * with optional markers that can be interactively moved.
 */
export default function Curve({ data, geometryKey }: { data: CurveData; geometryKey: string }) {
  const geometryDefaults = useAppStore((state) => state.geometryDefaults);

  // Merge with defaults from Pydantic (single source of truth)
  // Memoize to prevent recreating on every render
  const fullData = useMemo(
    () => getGeometryWithDefaults<CurveData>(data, "Curve", geometryDefaults),
    [data, geometryDefaults]
  );

  const {
    position: positionProp,
    color,
    material,
    thickness,
    divisions,
    marker,
    virtual_marker,
  } = fullData;

  const theme = useTheme();

  // Refs following Box.tsx pattern
  const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const hoverMeshRef = useRef<THREE.Mesh | null>(null);
  const drawingMarkerMeshRef = useRef<THREE.Mesh | null>(null);
  const virtualMarkerMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const curveRef = useRef<THREE.CatmullRomCurve3 | null>(null);

  const [markerCount, setMarkerCount] = useState(0);
  const [markerPositions, setMarkerPositions] = useState<THREE.Vector3[]>([]);
  const [lineSegments, setLineSegments] = useState<THREE.Vector3[]>([]);
  const [virtualMarkerPositions, setVirtualMarkerPositions] = useState<THREE.Vector3[]>([]);

  // Use individual selectors to prevent unnecessary re-renders
  const currentFrame = useAppStore((state) => state.currentFrame);
  const roomId = useAppStore((state) => state.roomId);
  const lock = useAppStore((state) => state.lock);
  const mode = useAppStore((state) => state.mode);
  const drawingIsValid = useAppStore((state) => state.drawingIsValid);
  const drawingPointerPosition = useAppStore((state) => state.drawingPointerPosition);
  const exitDrawingMode = useAppStore((state) => state.exitDrawingMode);
  const userName = useAppStore((state) => state.userName);
  const setGeometryFetching = useAppStore((state) => state.setGeometryFetching);
  const removeGeometryFetching = useAppStore((state) => state.removeGeometryFetching);
  const setCurveLength = useAppStore((state) => state.setCurveLength);
  const activeCurveForDrawing = useAppStore((state) => state.activeCurveForDrawing);
  const registerCurveRef = useAppStore((state) => state.registerCurveRef);
  const unregisterCurveRef = useAppStore((state) => state.unregisterCurveRef);
  const selections = useAppStore((state) => state.selections);
  const updateSelections = useAppStore((state) => state.updateSelections);
  const hoveredGeometryInstance = useAppStore((state) => state.hoveredGeometryInstance);
  const setHoveredGeometryInstance = useAppStore((state) => state.setHoveredGeometryInstance);
  const geometries = useAppStore((state) => state.geometries);
  const geometryUpdateSources = useAppStore((state) => state.geometryUpdateSources);

  const isDrawing = mode === 'drawing';
  const isEditing = mode === 'editing';

  // Fetch frame keys to check if required data is available
  const { data: frameKeysData, isLoading: isLoadingKeys } = useFrameKeys(roomId!, currentFrame);

  // Check if required keys are available for this geometry
  const requiredKeys = useMemo(() => {
    const keys: string[] = [];
    if (typeof positionProp === "string") keys.push(positionProp);
    return keys;
  }, [positionProp]);

  const hasRequiredKeys = useMemo(() => {
    // While loading keys, assume we have them (keep previous frame rendered)
    if (isLoadingKeys) return true;
    // If no keys data yet, assume we have them (keep previous frame)
    if (!frameKeysData?.keys) return true;
    // Only when we have keys data, check if required keys are available
    const availableKeys = new Set(frameKeysData.keys);
    return requiredKeys.every(key => availableKeys.has(key));
  }, [frameKeysData, requiredKeys, isLoadingKeys]);

  // Check if this curve is the active drawing target
  const isActiveDrawingTarget = activeCurveForDrawing === geometryKey;

  // Use geometry-specific selection (following Box.tsx pattern)
  const curveSelection = selections[geometryKey] || [];
  const selectionSet = useMemo(() => new Set(curveSelection), [curveSelection]);
  const selectedIndices = useMemo(() => Array.from(selectionSet), [selectionSet]);
  const validSelectedIndices = useMemo(
    () => selectedIndices.filter((id) => id < markerCount),
    [selectedIndices, markerCount]
  );

  // --- Data Fetching (following Box.tsx pattern) ---
  const { data: positionData, isFetching: isPositionFetching, isError: isPositionError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, positionProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [positionProp as string], signal),
    enabled: !!roomId && !!userName && typeof positionProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  // Check if any enabled query is still fetching
  const isFetching = typeof positionProp === "string" && isPositionFetching;

  // Check if any query has errored
  const hasQueryError = useMemo(
    () => typeof positionProp === "string" && isPositionError,
    [positionProp, isPositionError]
  );


  // Event handlers (following Box.tsx pattern)
  const handleMarkerClick = useCallback((event: any) => {
    if (event.detail !== 1 || event.instanceId === undefined) return;
    event.stopPropagation();
    updateSelections(geometryKey, event.instanceId, event.shiftKey);
  }, [geometryKey, updateSelections]);

  const onPointerEnterHandler = useCallback((event: any) => {
    if (event.instanceId === undefined) return;
    event.stopPropagation();
    setHoveredGeometryInstance(geometryKey, event.instanceId);
  }, [setHoveredGeometryInstance, geometryKey]);

  const onPointerLeaveHandler = useCallback(() => {
    setHoveredGeometryInstance(null, null);
  }, [setHoveredGeometryInstance]);

  // Check if this curve has static positions (can be edited)
  const hasStaticPosition = useMemo(() => {
    const currentGeometry = geometries[geometryKey];
    const position = currentGeometry?.data?.position;
    return Array.isArray(position) && position.length > 0 && Array.isArray(position[0]);
  }, [geometries, geometryKey]);

  const handleVirtualMarkerClick = useCallback((event: any) => {
    if (event.instanceId === undefined) return;
    event.stopPropagation();

    // Only allow virtual marker clicks in editing mode
    if (mode !== 'editing') return;

    const insertIndex = event.instanceId;
    const newPoints = [...markerPositions];
    newPoints.splice(insertIndex + 1, 0, virtualMarkerPositions[insertIndex].clone());

    // Convert to number[][] and update store
    const positions = newPoints.map(p => [p.x, p.y, p.z]);
    useAppStore.getState().updateGeometry(geometryKey, {
      type: "Curve",
      data: {
        ...fullData,
        position: positions,
      },
    }, 'local');

    // Select the newly added marker after refs update
    setTimeout(() => {
      updateSelections(geometryKey, insertIndex + 1, false);
    }, 0);
  }, [mode, markerPositions, virtualMarkerPositions, geometryKey, updateSelections, fullData]);

  const handleDrawingMarkerClick = useCallback(() => {
    if (drawingIsValid && drawingPointerPosition) {
      const newPoints = [...markerPositions];
      newPoints.push(drawingPointerPosition.clone());

      // Convert to number[][] and update store
      const positions = newPoints.map(p => [p.x, p.y, p.z]);
      useAppStore.getState().updateGeometry(geometryKey, {
        type: "Curve",
        data: {
          ...fullData,
          position: positions,
        },
      }, 'local');
    } else {
      exitDrawingMode();
    }
  }, [drawingIsValid, drawingPointerPosition, markerPositions, exitDrawingMode, geometryKey, fullData]);

  // Helper function to resolve color based on theme
  const getLineColor = useCallback((colorValue: string) => {
    return colorValue === "default"
      ? theme.palette.text.primary
      : colorValue;
  }, [theme]);

  const getMarkerColor = useCallback((colorValue: string) => {
    return colorValue === "default"
      ? theme.palette.secondary.main
      : colorValue;
  }, [theme]);

  // Report fetching state to global store
  useEffect(() => {
    setGeometryFetching(geometryKey, isFetching);
  }, [geometryKey, isFetching, setGeometryFetching]);

  // Clean up fetching state and curve ref on unmount
  useEffect(() => {
    return () => {
      removeGeometryFetching(geometryKey);
      unregisterCurveRef(geometryKey);
    };
  }, [geometryKey, removeGeometryFetching, unregisterCurveRef]);

  // Effect 1: Consolidated data processing and mesh update (following Box.tsx:265-390)
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

      // For curves with dynamic positions in editing mode:
      // If current data is unavailable but we have previous data, use previous positions
      // NOTE: markerPositions is intentionally NOT in dependencies to avoid infinite loops
      // We use the markerPositions value from closure which represents the last successful state
      const shouldPreserveState = isEditing && !hasStaticPosition && finalCount === 0 && markerCount > 0 && markerPositions.length === markerCount;

      if (finalCount === 0 && !shouldPreserveState) {
        // Normal case: reset to 0 if no data
        if (markerCount !== 0) setMarkerCount(0);
        return;
      }

      // Process position attribute using standardized utility
      // If preserving state, use stored markerPositions; otherwise use fetched data
      const finalPositions = shouldPreserveState
        ? markerPositions.flatMap(p => [p.x, p.y, p.z])
        : processPositionAttribute(positionProp, fetchedPosition);

      const actualCount = shouldPreserveState ? markerCount : finalCount;

      // --- Validation Step ---
      const isDataValid = validateArrayLengths(
        { positions: finalPositions },
        { positions: actualCount * 3 }
      );

      if (!isDataValid) {
        console.error("Curve data is invalid or has inconsistent lengths.");
        if (markerCount !== 0) setMarkerCount(0);
        return;
      }

      // --- Mesh Resizing Step ---
      if (markerCount !== actualCount) {
        setMarkerCount(actualCount);
        return;
      }

      // --- Main Marker Mesh Update ---
      const mainMesh = mainMeshRef.current;
      if (!mainMesh) return;

      const newMarkerPositions: THREE.Vector3[] = [];

      for (let i = 0; i < actualCount; i++) {
        const i3 = i * 3;
        _vec3.set(finalPositions[i3], finalPositions[i3 + 1], finalPositions[i3 + 2]);

        // Store for curve calculation
        newMarkerPositions.push(_vec3.clone());

        // Update instancedMesh (uniform scale for dodecahedron markers)
        _matrix.identity()
          .setPosition(_vec3)
          .scale(_vec3_2.set(marker.size, marker.size, marker.size));
        mainMesh.setMatrixAt(i, _matrix);
      }

      mainMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
      mainMesh.instanceMatrix.needsUpdate = true;

      mainMesh.computeBoundingBox();
      mainMesh.computeBoundingSphere();

      // --- Selection Mesh Update ---
      if (marker.selecting?.enabled && selectionMeshRef.current) {
        const selectionMesh = selectionMeshRef.current;
        validSelectedIndices.forEach((id, index) => {
          if (id >= actualCount) return;
          const i3 = id * 3;
          _vec3.set(finalPositions[i3], finalPositions[i3 + 1], finalPositions[i3 + 2]);

          // Apply SELECTION_SCALE
          const size = marker.size * SELECTION_SCALE;
          _matrix.identity()
            .setPosition(_vec3)
            .scale(_vec3_2.set(size, size, size));
          selectionMesh.setMatrixAt(index, _matrix);
        });
        selectionMesh.instanceMatrix.needsUpdate = true;
        selectionMesh.computeBoundingBox();
        selectionMesh.computeBoundingSphere();
      }

      // Store positions for curve calculation
      // When preserving state, keep existing markerPositions; otherwise update with new positions
      if (!shouldPreserveState) {
        // Only update if positions actually changed to avoid infinite loops
        const positionsChanged =
          markerPositions.length !== newMarkerPositions.length ||
          newMarkerPositions.some((pos, i) =>
            !markerPositions[i] ||
            !pos.equals(markerPositions[i])
          );

        if (positionsChanged) {
          setMarkerPositions(newMarkerPositions);
        }
      }

    } catch (error) {
      console.error("Error processing Curve data:", error);
      if (markerCount !== 0) setMarkerCount(0);
    }
  }, [
    data, // Add data to dependencies to ensure updates trigger
    isFetching,
    hasQueryError,
    positionData,
    positionProp,
    markerCount,
    validSelectedIndices,
    marker.selecting?.enabled,
    marker.size,
    geometryKey,
    isEditing, // Used in logic to preserve markerCount for dynamic positions
    hasStaticPosition, // Used in logic to preserve markerCount for dynamic positions
  ]);

  // Effect 2: Hover mesh updates (following Box.tsx:392-417)
  useEffect(() => {
    if (!marker.hovering?.enabled || !hoverMeshRef.current || !mainMeshRef.current) return;
    if (markerCount === 0) return;

    // Don't show hover in editing mode with dynamic positions
    if (isEditing && !hasStaticPosition) return;

    const hoverMesh = hoverMeshRef.current;
    const mainMesh = mainMeshRef.current;

    // Only show hover if it's for this geometry
    if (hoveredGeometryInstance?.geometryKey === geometryKey &&
        hoveredGeometryInstance?.instanceId !== null &&
        hoveredGeometryInstance.instanceId < markerCount) {
      hoverMesh.visible = true;

      // Get transform from main mesh
      mainMesh.getMatrixAt(hoveredGeometryInstance.instanceId, _matrix2);
      _matrix2.decompose(_vec3, _quat2, _vec3_2);

      // Apply HOVER_SCALE
      hoverMesh.position.copy(_vec3);
      hoverMesh.scale.copy(_vec3_2).multiplyScalar(HOVER_SCALE);
    } else {
      hoverMesh.visible = false;
    }
  }, [hoveredGeometryInstance, markerCount, marker.hovering?.enabled, geometryKey, isEditing, hasStaticPosition]);

  // Effect 3: Curve calculation from marker positions
  useEffect(() => {
    const allPoints: THREE.Vector3[] = [...markerPositions];
    if (!fullData) return;
    // Only add drawing pointer if this is the active drawing target
    if (drawingPointerPosition !== null && isDrawing && isActiveDrawingTarget) {
      allPoints.push(drawingPointerPosition);
    }
    if (allPoints.length < 2) {
      setLineSegments([]);
      setVirtualMarkerPositions([]);
      setCurveLength(0);
      // Unregister curve if we drop below 2 points
      curveRef.current = null;
      unregisterCurveRef(geometryKey);
      return;
    }
    const curve = new THREE.CatmullRomCurve3(allPoints, false, 'catmullrom', 0.5);
    curveRef.current = curve;

    // Register curve ref with store for Camera access
    registerCurveRef(geometryKey, curve);

    const curvePoints = curve.getPoints(divisions * allPoints.length);
    setLineSegments(curvePoints);

    setCurveLength(curve.getLength());

    // Only calculate virtual markers when they're actually used (editing mode with static positions)
    if (isEditing && hasStaticPosition && virtual_marker.enabled) {
      const _virtualMarkerPositions: THREE.Vector3[] = [];
      // virtual markers
      for (let i = 0; i < allPoints.length - 1; i++) {
        const controlPoint1 = allPoints[i];
        const controlPoint2 = allPoints[i + 1];

        // Find the curve points closest to each control point
        let minDist1 = Infinity;
        let minDist2 = Infinity;
        let index1 = 0;
        let index2 = 0;

        for (let j = 0; j < curvePoints.length; j++) {
          const dist1 = curvePoints[j].distanceTo(controlPoint1);
          const dist2 = curvePoints[j].distanceTo(controlPoint2);

          if (dist1 < minDist1) {
            minDist1 = dist1;
            index1 = j;
          }
          if (dist2 < minDist2) {
            minDist2 = dist2;
            index2 = j;
          }
        }

        // Place virtual marker at the midpoint index between these two curve points
        const midIndex = Math.floor((index1 + index2) / 2);
        _virtualMarkerPositions.push(curvePoints[midIndex].clone());
      }
      setVirtualMarkerPositions(_virtualMarkerPositions);
    } else {
      // Clear virtual markers when not in editing mode
      if (virtualMarkerPositions.length > 0) {
        setVirtualMarkerPositions([]);
      }
    }

  }, [markerPositions, drawingPointerPosition, mode, isActiveDrawingTarget, divisions, setCurveLength, geometryKey, registerCurveRef, unregisterCurveRef, isEditing, hasStaticPosition, virtual_marker.enabled, virtualMarkerPositions.length]);

  // Effect 4: Virtual markers mesh updates
  useEffect(() => {
    // Early returns - only log if we're actually doing work
    if (!isEditing || !hasStaticPosition || !virtual_marker.enabled) return;
    if (!virtualMarkerMeshRef.current) return;
    if (virtualMarkerPositions.length === 0) return;

    const virtualMesh = virtualMarkerMeshRef.current;

    for (let i = 0; i < virtualMarkerPositions.length; i++) {
      const pos = virtualMarkerPositions[i];
      _matrix.identity()
        .setPosition(pos)
        .scale(_vec3_2.set(virtual_marker.size, virtual_marker.size, virtual_marker.size));
      virtualMesh.setMatrixAt(i, _matrix);
    }

    virtualMesh.instanceMatrix.needsUpdate = true;
    virtualMesh.computeBoundingBox();
    virtualMesh.computeBoundingSphere();

  }, [virtualMarkerPositions, isEditing, hasStaticPosition, virtual_marker.enabled, virtual_marker.size]);

  // Effect 5: Drawing marker position and visibility update
  useEffect(() => {
    if (!drawingMarkerMeshRef.current) return;

    const drawingMesh = drawingMarkerMeshRef.current;

    // Set visibility based on drawing state
    if (!isDrawing || !isActiveDrawingTarget || !drawingPointerPosition || !marker.enabled) {
      drawingMesh.visible = false;
      // Reset scale and position when hiding to prevent ghost markers
      drawingMesh.scale.setScalar(0.001);
      drawingMesh.position.set(0, 0, 0);
      return;
    }

    drawingMesh.visible = true;
    drawingMesh.position.copy(drawingPointerPosition);

    // Set scale to match marker size
    drawingMesh.scale.setScalar(marker.size);

    // Color based on validity (red if invalid)
    const resolvedColor = marker.color === "default" ? theme.palette.secondary.main : marker.color;
    const materialColor = drawingIsValid ? resolvedColor : "#FF0000";
    (drawingMesh.material as THREE.MeshBasicMaterial).color.set(materialColor);

  }, [drawingPointerPosition, isDrawing, isActiveDrawingTarget, drawingIsValid, marker.color, marker.enabled, marker.size, theme.palette.secondary.main, geometryKey]);

  // Handle backspace/delete for selected curve markers in editing mode
  useEffect(() => {
    // Only handle deletion in editing mode
    if (mode !== 'editing' || selectedIndices.length === 0) {
      return;
    }

    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === "Backspace" || e.key === "Delete") {
        // Filter out selected indices
        const indicesToRemove = new Set(selectedIndices);
        const newPoints = markerPositions.filter((_, idx) => !indicesToRemove.has(idx));

        // Convert to number[][] and update store
        const positions = newPoints.map(p => [p.x, p.y, p.z]);
        useAppStore.getState().updateGeometry(geometryKey, {
          type: "Curve",
          data: {
            ...fullData,
            position: positions,
          },
        }, 'local');

        // Clear selection
        updateSelections(geometryKey, -1, false);
      }
    };

    window.addEventListener("keydown", handleKeyDown);
    return () => {
      window.removeEventListener("keydown", handleKeyDown);
    };
  }, [mode, selectedIndices, geometryKey, updateSelections, markerPositions, fullData]);

  // Persist position changes to server (debounced) - matches Box.tsx pattern
  const persistMarkerPositions = useCallback(async () => {
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
      await createGeometry(roomId, geometryKey, "Curve", currentGeometry.data, lock?.token);
    } catch (error: any) {
      console.error(`[Curve] Failed to persist ${geometryKey}:`, error);
    }
  }, [roomId, geometryKey, geometries, lock]);

  // Memoize debounced persist function to avoid recreation on every render
  const debouncedPersist = useMemo(
    () => debounce(persistMarkerPositions, 500),
    [persistMarkerPositions]
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

  // Get current position data for useGeometryEditing - simply convert markerPositions
  // markerPositions is already derived from store, so this is just a format conversion
  const finalPositionData = useMemo(() => {
    if (markerPositions.length === 0) return [];
    return markerPositions.map(point => [point.x, point.y, point.z]);
  }, [markerPositions]);

  // Use the geometry editing hook for transform controls
  useGeometryEditing(
    geometryKey,
    finalPositionData,
    selectedIndices,
    "Curve",
    fullData
  );

  // Shared geometry for all markers (following Box.tsx pattern)
  const markerGeometry = useMemo(() => {
    return new THREE.DodecahedronGeometry(1); // Unit size, scaled by matrix
  }, []);

  // Early returns after all hooks (Rules of Hooks)
  if (!roomId) return null;
  if (!marker || !virtual_marker) return null;
  // Don't render if geometry is disabled OR if required keys are not available
  if (fullData.active === false || !hasRequiredKeys) return null;

  // --- Render ---
  return (
    <group>
      {/* Curve line - UNIQUE TO CURVE */}
      {lineSegments.length >= 2 && (
        <Line
          points={lineSegments}
          color={isDrawing && isActiveDrawingTarget && !drawingIsValid ? "#FF0000" : getLineColor(color)}
          lineWidth={thickness}
          dashed={material === "LineDashedMaterial"}
        />
      )}

      {/* Main marker instancedMesh - following Box.tsx:494-506 */}
      {/* Always show markers except in specific cases where they should be hidden */}
      {marker.enabled && (
        <instancedMesh
          key={markerCount}
          ref={mainMeshRef}
          args={[undefined, undefined, markerCount]}
          // Only allow interaction if:
          // 1. Not in editing mode, OR
          // 2. In editing mode but position is static (can be transformed)
          onClick={marker.selecting?.enabled && (!isEditing || hasStaticPosition) ? handleMarkerClick : undefined}
          onPointerEnter={marker.hovering?.enabled && (!isEditing || hasStaticPosition) ? onPointerEnterHandler : undefined}
          onPointerLeave={marker.hovering?.enabled && (!isEditing || hasStaticPosition) ? onPointerLeaveHandler : undefined}
        >
          <primitive object={markerGeometry} attach="geometry" />
          <meshBasicMaterial
            color={getMarkerColor(marker.color)}
            opacity={marker.opacity}
            transparent={marker.opacity < 1}
          />
        </instancedMesh>
      )}

      {/* Selection mesh - following Box.tsx:509-523 */}
      {/* Only show selection if not in editing mode with dynamic positions */}
      {marker.enabled && marker.selecting?.enabled && (!isEditing || hasStaticPosition) && (
        <instancedMesh
          key={`selection-${validSelectedIndices.length}`}
          ref={selectionMeshRef}
          args={[undefined, undefined, validSelectedIndices.length]}
        >
          <primitive object={markerGeometry} attach="geometry" />
          <meshBasicMaterial
            side={THREE.FrontSide}
            transparent
            opacity={marker.selecting.opacity}
            color={marker.selecting.color}
          />
        </instancedMesh>
      )}

      {/* Hover mesh - following Box.tsx:526-536 */}
      {/* Only show hover if not in editing mode with dynamic positions */}
      {marker.enabled && marker.hovering?.enabled && (!isEditing || hasStaticPosition) && (
        <mesh ref={hoverMeshRef} visible={false}>
          <primitive object={markerGeometry} attach="geometry" />
          <meshBasicMaterial
            side={THREE.BackSide}
            transparent
            opacity={marker.hovering.opacity}
            color={marker.hovering.color}
          />
        </mesh>
      )}

      {/* Virtual markers instancedMesh - NEW (only in editing mode) */}
      {isEditing && hasStaticPosition && virtual_marker.enabled && (
        <instancedMesh
          key={`virtual-${virtualMarkerPositions.length}`}
          ref={virtualMarkerMeshRef}
          args={[undefined, undefined, virtualMarkerPositions.length]}
          onClick={handleVirtualMarkerClick}
        >
          <primitive object={markerGeometry} attach="geometry" />
          <meshBasicMaterial
            color={getMarkerColor(virtual_marker.color)}
            opacity={virtual_marker.opacity}
            transparent={virtual_marker.opacity < 1}
          />
        </instancedMesh>
      )}

      {/* Drawing marker - single mesh (only for active drawing target) */}
      {marker.enabled && isDrawing && isActiveDrawingTarget && drawingPointerPosition && (
        <mesh
          ref={drawingMarkerMeshRef}
          onClick={handleDrawingMarkerClick}
        >
          <primitive object={markerGeometry} attach="geometry" />
          <meshBasicMaterial
            color={getMarkerColor(marker.color)} // Color updated by effect
            opacity={marker.opacity}
            transparent={marker.opacity < 1}
          />
        </mesh>
      )}
    </group>
  );
}