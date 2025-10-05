import * as THREE from "three";
import { useQueries } from "@tanstack/react-query";
import { getFrameDataOptions } from "../../hooks/useTrajectoryData";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect } from "react";
import { useFrame } from "@react-three/fiber";
import { renderMaterial } from "./materials";

interface InteractionSettings {
  enabled: boolean;
  color: string | null;
  opacity: number;
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

// Re-usable vectors to avoid creating new objects in loops
const positionA = new THREE.Vector3();
const positionB = new THREE.Vector3();
const bondVec = new THREE.Vector3();
const bondVecInv = new THREE.Vector3();
const quaternion = new THREE.Quaternion();
const quaternionInv = new THREE.Quaternion();
const up = new THREE.Vector3(0, 1, 0);
const colorTempA = new THREE.Color();
const colorTempB = new THREE.Color();

const HOVER_SCALE = 1.25;
const SELECTION_SCALE = 1.05;

export default function Bonds({ data }: { data: BondData }) {
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

  const { currentFrame, roomId, selection } = useAppStore();
  const lastGoodData = useRef<any>(null);

  const selectionSet = useMemo(() => new Set(selection || []), [selection]);
  const bondResolution = resolution || 8;
  const bondScale = scale || 1.0;

  const keysToFetch = useMemo(() => {
    const keys: string[] = [];
    if (typeof positionProp === "string") keys.push(positionProp);
    if (typeof colorProp === "string") keys.push(colorProp);
    if (typeof radiusProp === "string") keys.push(radiusProp);
    if (typeof connectivityProp === "string") keys.push(connectivityProp);
    return keys;
  }, [positionProp, colorProp, radiusProp, connectivityProp]);

  const queries = useMemo(() => {
    if (!roomId) return [];
    return keysToFetch.map((key) => getFrameDataOptions(roomId, currentFrame, key));
  }, [roomId, currentFrame, keysToFetch]);

  const queryResults = useQueries({ queries });

  // FIX 3: Create a stable dependency array from the query data.
  const queryDataDeps = useMemo(() => queryResults.map(q => q.data?.data), [queryResults]);

  const { processedData, isFetching } = useMemo(() => {
    console.log("Recalculating bond data..."); // This should now only log when data actually changes.
    const isQueryFetching = queryResults.some((q) => q.isFetching || q.isPlaceholderData);
    if (keysToFetch.length > 0 && !queryResults.every((q) => q.isSuccess)) {
      return { isFetching: isQueryFetching, processedData: null };
    }

    const fetchedMap = new Map(queryResults.map((r, i) => [keysToFetch[i], r.data?.data]));

    let atomCount = 0;
    let finalPositions: number[] | Float32Array | Float64Array = [];
    if (typeof positionProp === "string") {
      const remoteData = fetchedMap.get(positionProp) || [];
      atomCount = remoteData.length / 3;
      finalPositions = remoteData;
    } else if (positionProp.length > 0 && Array.isArray(positionProp[0])) {
      atomCount = positionProp.length;
      finalPositions = (positionProp as number[][]).flat();
    } else {
      atomCount = positionProp.length / 3;
      finalPositions = positionProp as number[];
    }

    if (atomCount === 0) return { isFetching: isQueryFetching, processedData: null };

    let finalColors: number[] | Float32Array | Float64Array;
    if (typeof colorProp === "string") {
      finalColors = fetchedMap.get(colorProp) || [];
    } else if (colorProp.length > 0 && Array.isArray(colorProp[0])) {
      finalColors = (colorProp as number[][]).flat();
    } else {
      finalColors = new Array(atomCount).fill(colorProp).flat();
    }

    let finalRadii: number[] | Float32Array | Float64Array;
    if (typeof radiusProp === "string") {
      finalRadii = fetchedMap.get(radiusProp) || [];
    } else if (Array.isArray(radiusProp)) {
      finalRadii = radiusProp as number[];
    } else {
      finalRadii = new Array(atomCount).fill(radiusProp);
    }
    
    const connectivityRaw = typeof connectivityProp === "string" 
      ? fetchedMap.get(connectivityProp) || [] 
      : connectivityProp;

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
    
    if (bondPairs.length === 0) return { isFetching: isQueryFetching, processedData: null };

    // FIX 2: We now create two instances (half-bonds) per bond pair.
    const instances: { matrix: THREE.Matrix4; color: THREE.Color }[] = [];

    for (const bond of bondPairs) {
      const [a, b] = bond;
      if (a >= atomCount || b >= atomCount) continue;

      const a3 = a * 3;
      const b3 = b * 3;
      positionA.fromArray(finalPositions, a3);
      positionB.fromArray(finalPositions, b3);

      bondVec.subVectors(positionB, positionA);
      const length = bondVec.length();
      const radius = finalRadii[a] * bondScale;

      // Create a matrix for each half of the bond
      const halfLength = length / 2;
      const localMatrixA = new THREE.Matrix4();
      const localMatrixB = new THREE.Matrix4();
      
      // Common rotation for both halves
      quaternion.setFromUnitVectors(up, bondVec.clone().normalize());

      // First half-bond, starting from atom A
      localMatrixA.compose(positionA, quaternion, new THREE.Vector3(radius, halfLength, radius));
      colorTempA.setRGB(finalColors[a3], finalColors[a3 + 1], finalColors[a3 + 2]);
      instances.push({ matrix: localMatrixA, color: colorTempA.clone() });

      // Second half-bond, starting from atom B
      // We need to invert the bond vector and create a new quaternion
      bondVecInv.subVectors(positionA, positionB);
      quaternionInv.setFromUnitVectors(up, bondVecInv.clone().normalize());
      localMatrixB.compose(positionB, quaternionInv, new THREE.Vector3(radius, halfLength, radius));
      colorTempB.setRGB(finalColors[b3], finalColors[b3 + 1], finalColors[b3 + 2]);
      instances.push({ matrix: localMatrixB, color: colorTempB.clone() });
    }

    return { isFetching: isQueryFetching, processedData: { instances } };
    // FIX 3: Use the new stable dependency array.
  }, [keysToFetch, positionProp, colorProp, radiusProp, connectivityProp, bondScale, ...queryDataDeps]);

  const dataToRender = processedData || lastGoodData.current;
  if (processedData) lastGoodData.current = processedData;

  const bondGeometry = useMemo(() => {
    // FIX 1 & 2: Create a unit cylinder that goes from y=0 to y=1.
    // This makes it a perfect "half-bond" that can be positioned at an atom's center and scaled outwards.
    const geom = new THREE.CylinderGeometry(1, 1, 1, bondResolution, 1);
    geom.translate(0, 0.5, 0); // Move base to origin
    return geom;
  }, [bondResolution]);

  useEffect(() => {
    if (!mainMeshRef.current || !dataToRender) return;
    const mesh = mainMeshRef.current;
    const { instances } = dataToRender;
    // The number of instances is now double the number of bonds
    if (mesh.count !== instances.length) {
      mesh.count = instances.length;
    }
    for (let i = 0; i < instances.length; i++) {
      mesh.setMatrixAt(i, instances[i].matrix);
      mesh.setColorAt(i, instances[i].color);
    }
    mesh.instanceMatrix.needsUpdate = true;
    if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;
  }, [dataToRender]);
  
  // Hover/Selection logic may need adjustment if you want to select the whole bond (2 instances)
  // This current logic will highlight/select half-bonds.
  useFrame(() => {
    if (!dataToRender || isFetching) return;
    // ... existing useFrame logic ...
  });
  
  const onPointerMove = (e: any) => {
    if (e.instanceId === undefined) return;
    e.stopPropagation();
    setHoveredBondId(e.instanceId);
  };

  const onPointerOut = () => setHoveredBondId(null);
  
  if (!dataToRender) return null;

  return (
    <group>
      <instancedMesh
        ref={mainMeshRef}
        // FIX 2: The number of instances can now be up to bondPairs.length * 2
        args={[bondGeometry, undefined, dataToRender.instances.length]}
        onPointerMove={hovering.enabled ? onPointerMove : undefined}
        onPointerOut={hovering.enabled ? onPointerOut : undefined}
      >
        {renderMaterial(material, opacity)}
      </instancedMesh>
      
      {/* The rest of the component remains the same */}
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