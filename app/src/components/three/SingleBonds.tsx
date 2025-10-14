import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { getFrames } from "../../myapi/client";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect, useCallback } from "react";
import { renderMaterial } from "./materials";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import {
  processNumericAttribute,
  processColorAttribute,
  SELECTION_SCALE,
  HOVER_SCALE,
} from "../../utils/geometryData";
import { _vec3, _vec3_2, _quat, _matrix, _color } from "../../utils/threeObjectPools";
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
  color: string | number[][] | string;
  radius: string | number[] | number;
  material: string;
  resolution: number;
  scale: number;
  opacity: number;
  selecting: InteractionSettings;
  hovering: InteractionSettings;
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

  const { roomId, currentFrame, clientId, selections, updateSelections, setGeometryFetching, removeGeometryFetching, requestPathtracingUpdate } = useAppStore();

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

  const { data: connectivityData, isFetching: isConnectivityFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, connectivityProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [connectivityProp as string], signal),
    enabled: !!roomId && !!clientId && typeof connectivityProp === "string",
    placeholderData: keepPreviousData,
  });

  // Check if any enabled query is still fetching
  const isFetching =
    (typeof positionProp === "string" && isPositionFetching) ||
    (typeof colorProp === "string" && shouldFetchAsFrameData(colorProp as string) && isColorFetching) ||
    (typeof radiusProp === "string" && isRadiusFetching) ||
    (typeof connectivityProp === "string" && isConnectivityFetching);

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
      // Process atom positions first to get atom count
      const fetchedPosition = typeof positionProp === 'string' ? positionData?.[positionProp as string] : undefined;
      let atomCount = 0;
      let finalPositions: number[] = [];

      if (fetchedPosition) {
        atomCount = fetchedPosition.length / 3;
        finalPositions = Array.from(fetchedPosition);
      } else if (typeof positionProp !== "string") {
        if (Array.isArray(positionProp) && Array.isArray(positionProp[0])) {
          atomCount = (positionProp as number[][]).length;
          finalPositions = (positionProp as number[][]).flat();
        } else if (Array.isArray(positionProp)) {
          atomCount = (positionProp as number[]).length / 3;
          finalPositions = positionProp as number[];
        }
      }

      if (atomCount === 0) {
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      // Process colors
      const fetchedColor = typeof colorProp === 'string' ? colorData?.[colorProp as string] : undefined;
      const finalColors = processColorAttribute(colorProp, fetchedColor, atomCount);

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

        const bondVec = _vec3_2.clone().sub(_vec3);
        const length = bondVec.length();
        const radius = finalRadii[a] * bondScale;
        const halfLength = length / 2;

        // First half-bond (from atom A)
        _quat.setFromUnitVectors(_up, bondVec.clone().normalize());
        _matrix.compose(_vec3, _quat, new THREE.Vector3(radius, halfLength, radius));
        mainMesh.setMatrixAt(instanceIndex, _matrix);
        _colorA.setRGB(finalColors[a3], finalColors[a3 + 1], finalColors[a3 + 2]);
        mainMesh.setColorAt(instanceIndex, _colorA);
        instanceIndex++;

        // Second half-bond (from atom B, inverted direction)
        const bondVecInv = _vec3.clone().sub(_vec3_2);
        _quatInv.setFromUnitVectors(_up, bondVecInv.clone().normalize());
        _matrix.compose(_vec3_2, _quatInv, new THREE.Vector3(radius, halfLength, radius));
        mainMesh.setMatrixAt(instanceIndex, _matrix);
        _colorB.setRGB(finalColors[b3], finalColors[b3 + 1], finalColors[b3 + 2]);
        mainMesh.setColorAt(instanceIndex, _colorB);
        instanceIndex++;
      }

      mainMesh.instanceMatrix.needsUpdate = true;
      if (mainMesh.instanceColor) mainMesh.instanceColor.needsUpdate = true;

      // --- Selection Mesh Update ---
      if (selecting.enabled && selectionMeshRef.current) {
        const selectionMesh = selectionMeshRef.current;
        selectedBondIndices.forEach((bondInstanceId, arrayIndex) => {
          // Copy matrix from main mesh and scale up
          mainMesh.getMatrixAt(bondInstanceId, _matrix);
          _matrix.scale(new THREE.Vector3(SELECTION_SCALE, SELECTION_SCALE, SELECTION_SCALE));
          selectionMesh.setMatrixAt(arrayIndex, _matrix);
        });
        selectionMesh.instanceMatrix.needsUpdate = true;
      }

      // --- Hover Mesh Update ---
      if (hovering?.enabled && hoverMeshRef.current) {
        const hoverMesh = hoverMeshRef.current;
        if (hoveredBondId !== null && hoveredBondId < finalCount) {
          hoverMesh.visible = true;
          // Get the matrix from the hovered bond instance
          mainMesh.getMatrixAt(hoveredBondId, _matrix);
          // Decompose to get position, rotation, and scale
          const position = new THREE.Vector3();
          const rotation = new THREE.Quaternion();
          const scale = new THREE.Vector3();
          _matrix.decompose(position, rotation, scale);
          // Apply to hover mesh with scaled-up dimensions
          hoverMesh.position.copy(position);
          hoverMesh.quaternion.copy(rotation);
          hoverMesh.scale.copy(scale).multiplyScalar(HOVER_SCALE);
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