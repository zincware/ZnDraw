import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { getFrames, updateGeometryActive } from "../../myapi/client";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect, useCallback } from "react";
import { renderMaterial } from "./materials";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import {
  processNumericAttribute,
  processColorData,
  expandSharedColor,
  SELECTION_SCALE,
  HOVER_SCALE,
} from "../../utils/geometryData";
import { _vec3, _vec3_2, _vec3_3, _vec3_4, _quat, _quat2, _matrix, _matrix2, _color } from "../../utils/threeObjectPools";
import { convertInstancedMeshToMerged, disposeMesh } from "../../utils/convertInstancedMesh";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

interface InteractionSettings {
  enabled: boolean;
  color: string | null;
  opacity: number;
  is_drawable: boolean;
}

interface BondData {
  position: string | number[][];
  connectivity?: Array<[number, number, number]> | string;
  color: string | string[]; // Dynamic ref or list of hex strings
  radius: string | number[] | number;
  material: string;
  resolution: number;
  scale: number;
  opacity: number;
  selecting: InteractionSettings;
  hovering: InteractionSettings;
  active?: boolean; // Whether geometry is active (can be disabled on critical errors)
}

// Bond-specific reusable THREE objects
const _up = new THREE.Vector3(0, 1, 0);
const _colorA = new THREE.Color();
const _colorB = new THREE.Color();
const _quatInv = new THREE.Quaternion();

export default function Bonds({
  data,
  geometryKey,
  pathtracingEnabled = false
}: {
  data: BondData;
  geometryKey: string;
  pathtracingEnabled?: boolean;
}) {
  const { geometryDefaults } = useAppStore();

  // Merge with defaults from Pydantic (single source of truth)
  const fullData = getGeometryWithDefaults<BondData>(data, "Bond", geometryDefaults);

  const {
    position: positionProp,
    color: colorProp,
    radius: radiusProp,
    connectivity: connectivityProp,
    material,
    resolution,
    scale,
    opacity,
    selecting,
    hovering,
  } = fullData;

  const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const hoverMeshRef = useRef<THREE.Mesh | null>(null);
  const mergedMeshRef = useRef<THREE.Mesh | null>(null);
  const [hoveredBondId, setHoveredBondId] = useState<number | null>(null);
  const [instanceCount, setInstanceCount] = useState(0);
  const [bondPairs, setBondPairs] = useState<[number, number][]>([]);
  const hasDisabledGeometryRef = useRef(false);

  const { roomId, currentFrame, frameCount, clientId, selections, updateSelections, setGeometryFetching, removeGeometryFetching, requestPathtracingUpdate, updateGeometry, showSnackbar } = useAppStore();

  // Use geometry-specific selection
  const bondSelection = selections[geometryKey] || [];
  const selectionSet = useMemo(() => new Set(bondSelection), [bondSelection]);

  // Calculate which bond instances should be selected
  // Check if the bond pair index is in the selection set
  const selectedBondIndices = useMemo(() => {
    const indices: number[] = [];
    let instanceIndex = 0;
    for (let bondPairIndex = 0; bondPairIndex < bondPairs.length; bondPairIndex++) {
      const isSelected = selectionSet.has(bondPairIndex);
      if (isSelected) {
        indices.push(instanceIndex);     // First half-bond
        indices.push(instanceIndex + 1); // Second half-bond
      }
      instanceIndex += 2;
    }
    return indices;
  }, [bondPairs, selectionSet]);

  const bondResolution = resolution || 8;
  const bondScale = scale || 1.0;

  // Individual queries for each attribute - enables perfect cross-component caching
  const { data: positionData, isFetching: isPositionFetching, isError: isPositionError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, positionProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [positionProp as string], signal),
    enabled: !!roomId && !!clientId && frameCount > 0 && typeof positionProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: colorData, isFetching: isColorFetching, isError: isColorError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, colorProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [colorProp as string], signal),
    enabled: !!roomId && !!clientId && frameCount > 0 && typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string),
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: radiusData, isFetching: isRadiusFetching, isError: isRadiusError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, radiusProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [radiusProp as string], signal),
    enabled: !!roomId && !!clientId && frameCount > 0 && typeof radiusProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: connectivityData, isFetching: isConnectivityFetching, isError: isConnectivityError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, connectivityProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [connectivityProp as string], signal),
    enabled: !!roomId && !!clientId && frameCount > 0 && typeof connectivityProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  // Check if any enabled query is still fetching
  const isFetching =
    (typeof positionProp === "string" && isPositionFetching) ||
    (typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string) && isColorFetching) ||
    (typeof radiusProp === "string" && isRadiusFetching) ||
    (typeof connectivityProp === "string" && isConnectivityFetching);

  // Check if any query has errored - treat as data unavailable
  const hasQueryError = useMemo(
    () =>
      (typeof positionProp === "string" && isPositionError) ||
      (typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string) && isColorError) ||
      (typeof radiusProp === "string" && isRadiusError) ||
      (typeof connectivityProp === "string" && isConnectivityError),
    [positionProp, isPositionError, colorProp, isColorError, radiusProp, isRadiusError, connectivityProp, isConnectivityError]
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

  // Detect critical fetch failures and disable geometry
  useEffect(() => {
    if (!roomId || !clientId || hasDisabledGeometryRef.current || isFetching || frameCount === 0) {
      return;
    }

    // Critical errors for Bonds: position or connectivity queries failed
    const hasCriticalError =
      (typeof positionProp === "string" && isPositionError) ||
      (typeof connectivityProp === "string" && isConnectivityError);

    if (hasCriticalError) {
      if (process.env.NODE_ENV !== 'production') {
        console.warn(
          `Bonds geometry "${geometryKey}": Critical data fetch failed. Disabling geometry.`,
          {
            positionError: typeof positionProp === "string" && isPositionError,
            connectivityError: typeof connectivityProp === "string" && isConnectivityError
          }
        );
      }
      hasDisabledGeometryRef.current = true;

      // Optimistically update local state immediately
      const updatedGeometry = {
        type: "Bond",
        data: { ...data, active: false }
      };
      updateGeometry(geometryKey, updatedGeometry);

      // Show snackbar notification
      showSnackbar(`Geometry "${geometryKey}" disabled - data fetch failed`, "warning");

      // Then update server (server will skip emitting back to this client)
      updateGeometryActive(roomId, clientId, geometryKey, "Bond", false)
        .then(() => {
          if (process.env.NODE_ENV !== 'production') {
            console.info(`Bonds geometry "${geometryKey}" disabled successfully on server.`);
          }
        })
        .catch((error) => {
          console.error(`Failed to disable Bonds geometry "${geometryKey}" on server:`, error);
          // Rollback optimistic update on error
          const rollbackGeometry = {
            type: "Bond",
            data: { ...data }
          };
          updateGeometry(geometryKey, rollbackGeometry);
          hasDisabledGeometryRef.current = false;
        });
    }
  }, [
    roomId,
    clientId,
    geometryKey,
    frameCount,
    isFetching,
    positionProp,
    connectivityProp,
    isPositionError,
    isConnectivityError,
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
      // Process atom positions first to get atom count
      const fetchedPosition = typeof positionProp === 'string' ? positionData?.[positionProp as string] : undefined;

      // Calculate atom count from fetched or static position data
      let atomCount = 0;
      if (fetchedPosition) {
        atomCount = fetchedPosition.length / 3;
      } else if (Array.isArray(positionProp)) {
        atomCount = positionProp.length;  // positionProp is number[][]
      }

      if (atomCount === 0) {
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      // Use processNumericAttribute for consistent position processing
      const finalPositions = processNumericAttribute(positionProp, fetchedPosition, atomCount);

      // Process colors
      const fetchedColor = typeof colorProp === 'string' ? colorData?.[colorProp as string] : undefined;
      const colorHexArray = processColorData(colorProp, fetchedColor, atomCount);

      // Handle shared color (single color for all atoms)
      const finalColorHex = expandSharedColor(colorHexArray, atomCount);

      // Process radii
      const fetchedRadius = typeof radiusProp === 'string' ? radiusData?.[radiusProp as string] : undefined;
      const finalRadii = processNumericAttribute(radiusProp, fetchedRadius, atomCount);

      // Process connectivity
      const fetchedConnectivity = typeof connectivityProp === 'string' ? connectivityData?.[connectivityProp as string] : undefined;
      const connectivityRaw = fetchedConnectivity || connectivityProp;

      const newBondPairs: [number, number][] = [];
      if (connectivityRaw && connectivityRaw.length > 0) {
        if (Array.isArray(connectivityRaw[0])) {
          (connectivityRaw as number[][]).forEach(conn => newBondPairs.push([conn[0], conn[1]]));
        } else {
          for (let i = 0; i < connectivityRaw.length; i += 3) {
            newBondPairs.push([connectivityRaw[i], connectivityRaw[i + 1]]);
          }
        }
      }

      if (newBondPairs.length === 0) {
        if (instanceCount !== 0) setInstanceCount(0);
        if (bondPairs.length !== 0) setBondPairs([]);
        return;
      }

      // Each bond pair creates 2 half-bond instances
      const finalCount = newBondPairs.length * 2;

      // --- Mesh Resizing Step ---
      if (instanceCount !== finalCount) {
        setInstanceCount(finalCount);
        setBondPairs(newBondPairs);
        return;
      }

      // Update bondPairs if connectivity changed
      if (JSON.stringify(bondPairs) !== JSON.stringify(newBondPairs)) {
        setBondPairs(newBondPairs);
      }

      // --- Mesh Instance Update Step ---
      const mainMesh = mainMeshRef.current;
      if (!mainMesh) return;
      let instanceIndex = 0;
      for (const bond of newBondPairs) {
        const [a, b] = bond;
        if (a >= atomCount || b >= atomCount) continue;

        const a3 = a * 3;
        const b3 = b * 3;
        _vec3.fromArray(finalPositions, a3);
        _vec3_2.fromArray(finalPositions, b3);

        // Calculate bond vector using pooled object _vec3_3
        _vec3_3.copy(_vec3_2).sub(_vec3);
        const length = _vec3_3.length();
        const radius = finalRadii[a] * bondScale;
        const halfLength = length / 2;

        // Use _vec3_4 for scale vector
        _vec3_4.set(radius, halfLength, radius);

        // First half-bond (from atom A)
        // Normalize _vec3_3 in place instead of cloning
        _vec3_3.normalize();
        _quat.setFromUnitVectors(_up, _vec3_3);
        _matrix.compose(_vec3, _quat, _vec3_4);
        mainMesh.setMatrixAt(instanceIndex, _matrix);

        // Set color directly from hex string (THREE.Color.set() accepts hex)
        _colorA.set(finalColorHex[a]);
        mainMesh.setColorAt(instanceIndex, _colorA);
        instanceIndex++;

        // Second half-bond (from atom B, inverted direction)
        // Reuse _vec3_3 for inverted direction
        _vec3_3.copy(_vec3).sub(_vec3_2).normalize();
        _quatInv.setFromUnitVectors(_up, _vec3_3);
        _matrix.compose(_vec3_2, _quatInv, _vec3_4);
        mainMesh.setMatrixAt(instanceIndex, _matrix);

        // Set color directly from hex string (THREE.Color.set() accepts hex)
        _colorB.set(finalColorHex[b]);
        mainMesh.setColorAt(instanceIndex, _colorB);
        instanceIndex++;
      }

      mainMesh.instanceMatrix.needsUpdate = true;
      if (mainMesh.instanceColor) mainMesh.instanceColor.needsUpdate = true;

      // Update bounding box to prevent frustum culling issues
      mainMesh.computeBoundingBox();
      mainMesh.computeBoundingSphere();

      // --- Selection Mesh Update ---
      if (selecting.enabled && selectionMeshRef.current) {
        const selectionMesh = selectionMeshRef.current;
        selectedBondIndices.forEach((bondInstanceId, arrayIndex) => {
          // Copy matrix from main mesh and scale up using pooled vector
          mainMesh.getMatrixAt(bondInstanceId, _matrix);
          _vec3_3.set(SELECTION_SCALE, SELECTION_SCALE, SELECTION_SCALE);
          _matrix.scale(_vec3_3);
          selectionMesh.setMatrixAt(arrayIndex, _matrix);
        });
        selectionMesh.instanceMatrix.needsUpdate = true;

        // Update bounding box for selection mesh
        selectionMesh.computeBoundingBox();
        selectionMesh.computeBoundingSphere();
      }

      // --- Hover Mesh Update ---
      if (hovering?.enabled && hoverMeshRef.current) {
        const hoverMesh = hoverMeshRef.current;
        if (hoveredBondId !== null && hoveredBondId < finalCount) {
          hoverMesh.visible = true;
          // Get the matrix from the hovered bond instance using pooled objects
          mainMesh.getMatrixAt(hoveredBondId, _matrix2);
          _matrix2.decompose(_vec3, _quat2, _vec3_2);
          // Apply to hover mesh with scaled-up dimensions
          hoverMesh.position.copy(_vec3);
          hoverMesh.quaternion.copy(_quat2);
          hoverMesh.scale.copy(_vec3_2).multiplyScalar(HOVER_SCALE);
        } else {
          hoverMesh.visible = false;
        }
      }
    } catch (error) {
      console.error("Error processing Bonds data:", error);
      if (instanceCount !== 0) setInstanceCount(0);
    }
  }, [
    isFetching,
    hasQueryError,
    positionData,
    colorData,
    radiusData,
    connectivityData,
    positionProp,
    colorProp,
    radiusProp,
    connectivityProp,
    instanceCount,
    bondScale,
    selecting,
    selectedBondIndices,
    hovering,
    hoveredBondId,
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
    // DO NOT depend on positionData/colorData/radiusData/connectivityData/bondPairs/selections here!
    // That causes unnecessary rebuilds. instanceCount only changes AFTER mesh update completes.
  ]);

  const bondGeometry = useMemo(() => {
    // Create a unit cylinder that goes from y=0 to y=1.
    // This makes it a perfect "half-bond" that can be positioned at an atom's center and scaled outwards.
    const geom = new THREE.CylinderGeometry(1, 1, 1, bondResolution, 1);
    geom.translate(0, 0.5, 0); // Move base to origin
    return geom;
  }, [bondResolution]);

  const onClickHandler = useCallback((event: any) => {
    if (event.detail !== 1 || event.instanceId === undefined) return;
    event.stopPropagation();
    // Convert instance ID to bond pair index (each pair has 2 instances)
    const bondPairIndex = Math.floor(event.instanceId / 2);
    updateSelections(geometryKey, bondPairIndex, event.shiftKey);
  }, [updateSelections, geometryKey]);

  const onPointerEnter = useCallback((e: any) => {
    if (e.instanceId === undefined) return;
    e.stopPropagation();
    setHoveredBondId(e.instanceId);
  }, []);

  const onPointerMove = useCallback((e: any) => {
    if (e.instanceId === undefined) return;
    e.stopPropagation();
    setHoveredBondId(e.instanceId);
  }, []);

  const onPointerOut = useCallback(() => setHoveredBondId(null), []);

  if (!clientId || !roomId) return null;

  // Don't render if geometry is disabled
  if (fullData.active === false) {
    return null;
  }

  return (
    <group>
      {/* Instanced mesh - visible when NOT pathtracing */}
      {/* NOTE: Interactions (click, hover) disabled when pathtracing enabled */}
      <instancedMesh
        key={instanceCount}
        ref={mainMeshRef}
        args={[bondGeometry, undefined, instanceCount]}
        visible={!pathtracingEnabled}
        onClick={!pathtracingEnabled && selecting.enabled ? onClickHandler : undefined}
        onPointerEnter={!pathtracingEnabled && hovering.enabled ? onPointerEnter : undefined}
        onPointerMove={!pathtracingEnabled && hovering.enabled ? onPointerMove : undefined}
        onPointerOut={!pathtracingEnabled && hovering.enabled ? onPointerOut : undefined}
      >
        {renderMaterial(material, opacity)}
      </instancedMesh>

      {/* Selection and hover meshes - only when NOT pathtracing */}
      {!pathtracingEnabled && (
        <>
          {/* Selection mesh */}
          {selecting.enabled && (
            <instancedMesh
              key={`selection-${selectedBondIndices.length}`}
              ref={selectionMeshRef}
              args={[bondGeometry, undefined, selectedBondIndices.length]}
            >
              <meshBasicMaterial
                side={THREE.FrontSide}
                transparent
                opacity={selecting.opacity}
                color={selecting.color || "#FFFF00"}
              />
            </instancedMesh>
          )}

          {/* Hover mesh */}
          {hovering?.enabled && (
            <mesh ref={hoverMeshRef} visible={false}>
              <primitive object={bondGeometry} attach="geometry" />
              <meshBasicMaterial
                side={THREE.BackSide}
                transparent
                opacity={hovering.opacity}
                color={hovering.color || "#00FFFF"}
              />
            </mesh>
          )}
        </>
      )}

      {/* Merged mesh - visible when pathtracing */}
      {pathtracingEnabled && mergedMeshRef.current && (
        <primitive object={mergedMeshRef.current} />
      )}
    </group>
  );
}