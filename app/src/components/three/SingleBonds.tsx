import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { getFrames } from "../../myapi/client";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect } from "react";
import { renderMaterial } from "./materials";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import {
  type DataProp,
  processNumericAttribute,
  processColorAttribute,
  SELECTION_SCALE,
  HOVER_SCALE,
} from "../../utils/geometryData";
import { _vec3, _vec3_2, _quat, _matrix, _color } from "../../utils/threeObjectPools";

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

export default function Bonds({ data, geometryKey }: { data: BondData; geometryKey: string }) {
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
  } = data;

  const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const hoverMeshRef = useRef<THREE.Mesh | null>(null);
  const [hoveredBondId, setHoveredBondId] = useState<number | null>(null);
  const [instanceCount, setInstanceCount] = useState(0);

  const { roomId, currentFrame, clientId, selection, setGeometryFetching, removeGeometryFetching } = useAppStore();

  const selectionSet = useMemo(() => new Set(selection || []), [selection]);
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

      const bondPairs: [number, number][] = [];
      if (connectivityRaw && connectivityRaw.length > 0) {
        if (Array.isArray(connectivityRaw[0])) {
          (connectivityRaw as number[][]).forEach(conn => bondPairs.push([conn[0], conn[1]]));
        } else {
          for (let i = 0; i < connectivityRaw.length; i += 3) {
            bondPairs.push([connectivityRaw[i], connectivityRaw[i + 1]]);
          }
        }
      }

      if (bondPairs.length === 0) {
        if (instanceCount !== 0) setInstanceCount(0);
        return;
      }

      // Each bond pair creates 2 half-bond instances
      const finalCount = bondPairs.length * 2;

      // --- Mesh Resizing Step ---
      if (instanceCount !== finalCount) {
        setInstanceCount(finalCount);
        return;
      }

      // --- Mesh Instance Update Step ---
      const mainMesh = mainMeshRef.current;
      if (!mainMesh) return;

      let instanceIndex = 0;
      for (const bond of bondPairs) {
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
  ]);

  const bondGeometry = useMemo(() => {
    // Create a unit cylinder that goes from y=0 to y=1.
    // This makes it a perfect "half-bond" that can be positioned at an atom's center and scaled outwards.
    const geom = new THREE.CylinderGeometry(1, 1, 1, bondResolution, 1);
    geom.translate(0, 0.5, 0); // Move base to origin
    return geom;
  }, [bondResolution]);
  
  const onPointerMove = (e: any) => {
    if (e.instanceId === undefined) return;
    e.stopPropagation();
    setHoveredBondId(e.instanceId);
  };

  const onPointerOut = () => setHoveredBondId(null);

  if (!clientId || !roomId) return null;

  return (
    <group>
      <instancedMesh
        key={instanceCount}
        ref={mainMeshRef}
        args={[bondGeometry, undefined, instanceCount]}
        onPointerMove={hovering.enabled ? onPointerMove : undefined}
        onPointerOut={hovering.enabled ? onPointerOut : undefined}
      >
        {renderMaterial(material, opacity)}
      </instancedMesh>

      {/* Selection and hover meshes can be enabled later if needed */}
      {/* {selecting.enabled && (
        <instancedMesh ref={selectionMeshRef} args={[bondGeometry, undefined, 0]}>
          <meshBasicMaterial transparent opacity={selecting.opacity} color={selecting.color || "#FFFF00"} />
        </instancedMesh>
      )}
      {hovering.enabled && (
        <mesh ref={hoverMeshRef} visible={false} geometry={bondGeometry}>
          <meshBasicMaterial transparent opacity={hovering.opacity} color={hovering.color || "#00FFFF"} />
        </mesh>
      )} */}
    </group>
  );
}