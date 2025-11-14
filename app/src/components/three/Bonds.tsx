import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { getFrames } from "../../myapi/client";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect, useCallback, useReducer } from "react";
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
import { useFrameKeys } from "../../hooks/useSchemas";

interface InteractionSettings {
  enabled: boolean;
  color: string | null;
  opacity: number;
  is_drawable: boolean;
}

interface BondData {
  position: string | number[][];
  connectivity?: Array<[number, number, number | null]> | string;
  color: string | string[]; // Dynamic ref or list of hex strings
  radius: string | number[] | number;
  material: string;
  resolution: number;
  scale: number;
  opacity: number;
  bond_order_mode: "parallel" | "ignore";
  bond_order_offset: number;
  bond_order_radius_scale: { [key: string]: number };
  selecting: InteractionSettings;
  hovering: InteractionSettings;
  active?: boolean; // Whether geometry is active (can be disabled on critical errors)
}

// Bond-specific reusable THREE objects
const _up = new THREE.Vector3(0, 1, 0);
const _colorA = new THREE.Color();
const _colorB = new THREE.Color();
const _quatInv = new THREE.Quaternion();
const _perpVector = new THREE.Vector3();  // For perpendicular offsets
const _offsetPos = new THREE.Vector3();   // For offset positions
const _negatedDirection = new THREE.Vector3(); // For negated bond direction (avoid clone)
const _scaleVector = new THREE.Vector3(); // Reusable scale vector

/**
 * Calculate number of cylinders to render for a given bond order.
 * @param bondOrder The bond order (1, 1.5, 2, 3, or null)
 * @param mode The visualization mode ("parallel" or "ignore")
 * @returns Number of cylinders to render
 */
function getCylinderCount(bondOrder: number | null, mode: "parallel" | "ignore"): number {
  if (mode === "ignore") return 1;

  const order = bondOrder ?? 1;
  if (order === 1.5) return 2;  // Aromatic: 2 cylinders
  return Math.max(1, Math.round(order));  // 1, 2, or 3 cylinders
}

/**
 * Create a fast hash of bond pairs for comparison.
 * Much faster than JSON.stringify for large arrays.
 */
function hashBondPairs(pairs: Array<[number, number, number]>): string {
  let hash = '';
  for (let i = 0; i < pairs.length; i++) {
    const [a, b, o] = pairs[i];
    hash += `${a},${b},${o};`;
  }
  return hash;
}

// State for bond rendering
interface BondState {
  instanceCount: number;
  bondPairs: Array<[number, number, number]>;
  instanceToBondMap: number[];
  bondPairsHash: string;
}

type BondStateAction =
  | { type: 'UPDATE_ALL'; payload: Omit<BondState, 'bondPairsHash'> }
  | { type: 'RESET' };

function bondStateReducer(state: BondState, action: BondStateAction): BondState {
  switch (action.type) {
    case 'UPDATE_ALL':
      return {
        ...action.payload,
        bondPairsHash: hashBondPairs(action.payload.bondPairs),
      };
    case 'RESET':
      return { instanceCount: 0, bondPairs: [], instanceToBondMap: [], bondPairsHash: '' };
    default:
      return state;
  }
}

export default function Bonds({
  data,
  geometryKey,
  pathtracingEnabled = false
}: {
  data: BondData;
  geometryKey: string;
  pathtracingEnabled?: boolean;
}) {
  const geometryDefaults = useAppStore((state) => state.geometryDefaults);

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
    bond_order_mode,
    bond_order_offset,
    bond_order_radius_scale,
    selecting,
    hovering,
  } = fullData;

  const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const hoverMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const mergedMeshRef = useRef<THREE.Mesh | null>(null);
  const [hoveredBondPairIndex, setHoveredBondPairIndex] = useState<number | null>(null);

  // Consolidated state with useReducer (better performance for related state updates)
  const [bondState, dispatchBondState] = useReducer(bondStateReducer, {
    instanceCount: 0,
    bondPairs: [],
    instanceToBondMap: [],
    bondPairsHash: '',
  });
  const { instanceCount, bondPairs, instanceToBondMap, bondPairsHash } = bondState;

  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const currentFrame = useAppStore((state) => state.currentFrame);
  const frameCount = useAppStore((state) => state.frameCount);
  const userName = useAppStore((state) => state.userName);
  const selections = useAppStore((state) => state.selections);
  const updateSelections = useAppStore((state) => state.updateSelections);
  const setGeometryFetching = useAppStore((state) => state.setGeometryFetching);
  const removeGeometryFetching = useAppStore((state) => state.removeGeometryFetching);
  const requestPathtracingUpdate = useAppStore((state) => state.requestPathtracingUpdate);

  // Fetch frame keys to check if required data is available
  const { data: frameKeysData, isLoading: isLoadingKeys } = useFrameKeys(roomId!, currentFrame);

  // Check if required keys are available for this geometry
  const requiredKeys = useMemo(() => {
    const keys: string[] = [];
    if (typeof positionProp === "string") keys.push(positionProp);
    if (typeof connectivityProp === "string") keys.push(connectivityProp);
    return keys;
  }, [positionProp, connectivityProp]);

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
  const bondSelection = selections[geometryKey] || [];
  const selectionSet = useMemo(() => new Set(bondSelection), [bondSelection]);

  // Pre-compute instance ranges for each bond pair (memoized)
  const bondInstanceRanges = useMemo(() => {
    const ranges: Array<{ start: number; count: number }> = [];
    let instanceIndex = 0;
    for (let bondPairIndex = 0; bondPairIndex < bondPairs.length; bondPairIndex++) {
      const [, , bondOrder] = bondPairs[bondPairIndex];
      const cylinderCount = getCylinderCount(bondOrder, bond_order_mode);
      const instancesForThisBond = cylinderCount * 2;  // 2 half-bonds per cylinder
      ranges.push({ start: instanceIndex, count: instancesForThisBond });
      instanceIndex += instancesForThisBond;
    }
    return ranges;
  }, [bondPairs, bond_order_mode]);

  // Calculate which bond instances should be selected (optimized to only iterate selected bonds)
  const selectedBondIndices = useMemo(() => {
    if (bondSelection.length === 0) return [];

    const indices: number[] = [];
    // Only iterate through selected bonds, not all bonds
    for (const bondPairIndex of bondSelection) {
      const range = bondInstanceRanges[bondPairIndex];
      if (!range) continue; // Safety check
      // Add all instances for this selected bond
      for (let i = 0; i < range.count; i++) {
        indices.push(range.start + i);
      }
    }
    return indices;
  }, [bondSelection, bondInstanceRanges]);

  // Calculate which bond instances should be hovered (optimized using pre-computed ranges)
  const hoveredBondIndices = useMemo(() => {
    if (hoveredBondPairIndex === null) return [];

    const range = bondInstanceRanges[hoveredBondPairIndex];
    if (!range) return []; // Safety check

    const indices: number[] = [];
    for (let i = 0; i < range.count; i++) {
      indices.push(range.start + i);
    }
    return indices;
  }, [hoveredBondPairIndex, bondInstanceRanges]);

  const bondResolution = resolution || 8;
  const bondScale = scale || 1.0;

  // Individual queries for each attribute - enables perfect cross-component caching
  const { data: positionData, isFetching: isPositionFetching, isError: isPositionError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, positionProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [positionProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof positionProp === "string",
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

  const { data: radiusData, isFetching: isRadiusFetching, isError: isRadiusError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, radiusProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [radiusProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof radiusProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
  });

  const { data: connectivityData, isFetching: isConnectivityFetching, isError: isConnectivityError } = useQuery({
    queryKey: ["frame", roomId, currentFrame, connectivityProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [connectivityProp as string], signal),
    enabled: !!roomId && !!userName && frameCount > 0 && typeof connectivityProp === "string",
    placeholderData: keepPreviousData,
    retry: false,
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
        if (instanceCount !== 0) {
          dispatchBondState({ type: 'RESET' });
        }
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

      // Process connectivity with bond order
      const fetchedConnectivity = typeof connectivityProp === 'string' ? connectivityData?.[connectivityProp as string] : undefined;
      const connectivityRaw = fetchedConnectivity || connectivityProp;

      const newBondPairs: Array<[number, number, number]> = [];
      if (connectivityRaw && connectivityRaw.length > 0) {
        if (Array.isArray(connectivityRaw[0])) {
          // Array format: [[atom_a, atom_b, bond_order], ...]
          (connectivityRaw as Array<[number, number, number | null]>).forEach(conn => {
            const bondOrder = conn[2] ?? 1;  // Default to 1 if null
            newBondPairs.push([conn[0], conn[1], bondOrder]);
          });
        } else {
          // Flat array format: [atom_a, atom_b, bond_order, ...]
          // connectivityRaw might be a TypedArray or regular array
          for (let i = 0; i < connectivityRaw.length; i += 3) {
            const atomA = Number(connectivityRaw[i]);
            const atomB = Number(connectivityRaw[i + 1]);
            const bondOrder = Number(connectivityRaw[i + 2] ?? 1);  // Default to 1 if null
            newBondPairs.push([atomA, atomB, bondOrder]);
          }
        }
      }

      if (newBondPairs.length === 0) {
        if (instanceCount !== 0) {
          dispatchBondState({ type: 'RESET' });
        }
        return;
      }

      // Calculate total instance count based on bond orders and mode
      // Each bond pair creates N cylinders, each cylinder has 2 half-bonds
      // Also build mapping from instance index to bond pair index
      let totalInstances = 0;
      const newInstanceToBondMap: number[] = [];
      for (let bondPairIndex = 0; bondPairIndex < newBondPairs.length; bondPairIndex++) {
        const [, , bondOrder] = newBondPairs[bondPairIndex];
        const cylinderCount = getCylinderCount(bondOrder, bond_order_mode);
        const instancesForThisBond = cylinderCount * 2;  // 2 half-bonds per cylinder

        // Map all instances of this bond to the same bond pair index
        for (let i = 0; i < instancesForThisBond; i++) {
          newInstanceToBondMap.push(bondPairIndex);
        }

        totalInstances += instancesForThisBond;
      }
      const finalCount = totalInstances;

      // --- Mesh Resizing Step ---
      // Calculate hash for new bond pairs
      const newBondPairsHash = hashBondPairs(newBondPairs);

      // Only dispatch if values differ from current state
      // This prevents unnecessary updates while avoiding infinite loops
      // (instanceCount and bondPairsHash are NOT in deps, so state updates don't re-trigger this effect)
      if (finalCount !== instanceCount || newBondPairsHash !== bondPairsHash) {
        dispatchBondState({
          type: 'UPDATE_ALL',
          payload: {
            instanceCount: finalCount,
            bondPairs: newBondPairs,
            instanceToBondMap: newInstanceToBondMap,
          },
        });
        return;
      }

      // --- Mesh Instance Update Step ---
      const mainMesh = mainMeshRef.current;
      if (!mainMesh) return;
      let instanceIndex = 0;

      for (const bond of newBondPairs) {
        const [a, b, bondOrder] = bond;
        if (a >= atomCount || b >= atomCount) continue;

        const a3 = a * 3;
        const b3 = b * 3;
        _vec3.fromArray(finalPositions, a3);
        _vec3_2.fromArray(finalPositions, b3);

        // Calculate bond vector and normalize
        _vec3_3.copy(_vec3_2).sub(_vec3);
        const length = _vec3_3.length();
        _vec3_3.normalize();  // Bond direction vector

        // Get number of cylinders for this bond
        const cylinderCount = getCylinderCount(bondOrder, bond_order_mode);

        // Get radius scale for this bond order
        // Python dict keys like {2.0: 0.1} serialize to {"2.0": 0.1} in JSON
        // Try both "2" and "2.0" formats
        const bondOrderKey1 = String(bondOrder);
        const bondOrderKey2 = Number(bondOrder).toFixed(1);
        const radiusScale = bond_order_mode === "parallel"
          ? (bond_order_radius_scale[bondOrderKey1] ?? bond_order_radius_scale[bondOrderKey2] ?? 1.0)
          : 1.0;
        const radius = finalRadii[a] * bondScale * radiusScale;
        const halfLength = length / 2;

        // Pre-calculate perpendicular vector ONCE per bond (not per cylinder)
        // Only needed if multiple cylinders
        let perpVectorCalculated = false;
        if (cylinderCount > 1) {
          // Find a perpendicular vector to the bond direction
          // Use cross product with up vector, or right vector if bond is vertical
          if (Math.abs(_vec3_3.dot(_up)) < 0.99) {
            _perpVector.crossVectors(_vec3_3, _up).normalize();
          } else {
            // Bond is nearly vertical, use right vector instead
            _perpVector.set(1, 0, 0).cross(_vec3_3).normalize();
          }
          perpVectorCalculated = true;
        }

        // Render multiple cylinders for double/triple bonds
        for (let cylIdx = 0; cylIdx < cylinderCount; cylIdx++) {
          // Calculate offset for this cylinder
          // bond_order_offset is the spacing between adjacent cylinders
          let offsetAmount = 0;
          if (cylinderCount === 2) {
            // Double bond: cylinders at -0.5 and +0.5 spacing
            offsetAmount = (cylIdx === 0 ? -0.5 : 0.5) * bond_order_offset * radius;
          } else if (cylinderCount === 3) {
            // Triple bond: cylinders at -1, 0, +1 spacing (distance A-B = B-C = bond_order_offset * radius)
            offsetAmount = (cylIdx - 1) * bond_order_offset * radius;
          }

          // First half-bond (from atom A)
          if (perpVectorCalculated) {
            _offsetPos.copy(_vec3).addScaledVector(_perpVector, offsetAmount);
          } else {
            _offsetPos.copy(_vec3);
          }

          _vec3_4.set(radius, halfLength, radius);
          _quat.setFromUnitVectors(_up, _vec3_3);
          _matrix.compose(_offsetPos, _quat, _vec3_4);
          mainMesh.setMatrixAt(instanceIndex, _matrix);
          _colorA.set(finalColorHex[a]);
          mainMesh.setColorAt(instanceIndex, _colorA);
          instanceIndex++;

          // Second half-bond (from atom B, inverted direction)
          if (perpVectorCalculated) {
            _offsetPos.copy(_vec3_2).addScaledVector(_perpVector, offsetAmount);
          } else {
            _offsetPos.copy(_vec3_2);
          }

          // Inverted direction for second half (use pooled vector instead of clone)
          _vec3_4.set(radius, halfLength, radius);
          _negatedDirection.copy(_vec3_3).negate();
          _quatInv.setFromUnitVectors(_up, _negatedDirection);
          _matrix.compose(_offsetPos, _quatInv, _vec3_4);
          mainMesh.setMatrixAt(instanceIndex, _matrix);
          _colorB.set(finalColorHex[b]);
          mainMesh.setColorAt(instanceIndex, _colorB);
          instanceIndex++;
        }
      }

      mainMesh.instanceMatrix.needsUpdate = true;
      if (mainMesh.instanceColor) mainMesh.instanceColor.needsUpdate = true;

      // Update bounding box to prevent frustum culling issues
      mainMesh.computeBoundingBox();
      mainMesh.computeBoundingSphere();
    } catch (error) {
      console.error("Error processing Bonds data:", error);
      if (instanceCount !== 0) {
        dispatchBondState({ type: 'RESET' });
      }
    }
  }, [
    isFetching,
    // Data dependencies - these change when frame changes
    positionData,
    colorData,
    radiusData,
    connectivityData,
    // Prop dependencies
    positionProp,
    colorProp,
    radiusProp,
    connectivityProp,
    // Configuration
    bondScale,
    bond_order_mode,
    bond_order_offset,
    bond_order_radius_scale,
    // Note: instanceCount and bondPairsHash are NOT in dependencies to prevent infinite loops.
    // They are derived state - we compare calculated values against current state values
    // and only dispatch if different. Since they're not in deps, state updates don't
    // re-trigger this effect. This breaks the feedback loop while maintaining proper
    // React data flow (no refs needed).
  ]);

  // Separate effect for selection mesh updates (doesn't reprocess main mesh data)
  useEffect(() => {
    if (!selecting.enabled || !mainMeshRef.current || !selectionMeshRef.current) return;
    if (selectedBondIndices.length === 0) return;

    const mainMesh = mainMeshRef.current;
    const selectionMesh = selectionMeshRef.current;

    // Create scale vector once outside loop
    _scaleVector.set(SELECTION_SCALE, SELECTION_SCALE, SELECTION_SCALE);
    for (let arrayIndex = 0; arrayIndex < selectedBondIndices.length; arrayIndex++) {
      const bondInstanceId = selectedBondIndices[arrayIndex];
      mainMesh.getMatrixAt(bondInstanceId, _matrix);
      _matrix.scale(_scaleVector);
      selectionMesh.setMatrixAt(arrayIndex, _matrix);
    }
    selectionMesh.instanceMatrix.needsUpdate = true;

    // Only recompute bounding box if selection mesh instance count changed
    if (selectionMesh.count !== selectedBondIndices.length) {
      selectionMesh.computeBoundingBox();
      selectionMesh.computeBoundingSphere();
    }
  }, [selecting.enabled, selectedBondIndices]);

  // Separate effect for hover mesh updates (doesn't reprocess main mesh data)
  useEffect(() => {
    if (!hovering?.enabled || !mainMeshRef.current || !hoverMeshRef.current) return;
    if (hoveredBondIndices.length === 0) return;

    const mainMesh = mainMeshRef.current;
    const hoverMesh = hoverMeshRef.current;

    // Create scale vector once outside loop
    _scaleVector.set(HOVER_SCALE, HOVER_SCALE, HOVER_SCALE);
    for (let arrayIndex = 0; arrayIndex < hoveredBondIndices.length; arrayIndex++) {
      const bondInstanceId = hoveredBondIndices[arrayIndex];
      mainMesh.getMatrixAt(bondInstanceId, _matrix);
      _matrix.scale(_scaleVector);
      hoverMesh.setMatrixAt(arrayIndex, _matrix);
    }
    hoverMesh.instanceMatrix.needsUpdate = true;

    // Only recompute bounding box if hover mesh instance count changed
    if (hoverMesh.count !== hoveredBondIndices.length) {
      hoverMesh.computeBoundingBox();
      hoverMesh.computeBoundingSphere();
    }
  }, [hovering?.enabled, hoveredBondIndices]);

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
    // Convert instance ID to bond pair index using mapping
    // (handles bonds with different numbers of cylinders)
    const bondPairIndex = instanceToBondMap[event.instanceId] ?? Math.floor(event.instanceId / 2);
    updateSelections(geometryKey, bondPairIndex, event.shiftKey);
  }, [updateSelections, geometryKey, instanceToBondMap]);

  const onPointerEnter = useCallback((e: any) => {
    if (e.instanceId === undefined) return;
    e.stopPropagation();
    // Convert instance ID to bond pair index
    const bondPairIndex = instanceToBondMap[e.instanceId] ?? Math.floor(e.instanceId / 2);
    setHoveredBondPairIndex(bondPairIndex);
  }, [instanceToBondMap]);

  const onPointerMove = useCallback((e: any) => {
    if (e.instanceId === undefined) return;
    e.stopPropagation();
    // Convert instance ID to bond pair index
    const bondPairIndex = instanceToBondMap[e.instanceId] ?? Math.floor(e.instanceId / 2);
    setHoveredBondPairIndex(bondPairIndex);
  }, [instanceToBondMap]);

  const onPointerOut = useCallback(() => setHoveredBondPairIndex(null), []);

  if (!userName || !roomId) return null;

  // Don't render if geometry is disabled OR if required keys are not available
  if (fullData.active === false || !hasRequiredKeys) {
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
            <instancedMesh
              key={`hover-${hoveredBondIndices.length}`}
              ref={hoverMeshRef}
              args={[bondGeometry, undefined, hoveredBondIndices.length]}
            >
              <meshBasicMaterial
                side={THREE.BackSide}
                transparent
                opacity={hovering.opacity}
                color={hovering.color || "#00FFFF"}
              />
            </instancedMesh>
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