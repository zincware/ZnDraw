import * as THREE from "three";
import { useQueries } from "@tanstack/react-query";
import { getFrameDataOptions } from "../../hooks/useTrajectoryData";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState, useEffect, useCallback } from "react";
import { renderMaterial } from "./materials";
import { shouldFetchAsFrameData, hexToRgb } from "../../utils/colorUtils";

interface InteractionSettings {
  enabled: boolean;
  color: string;
  opacity: number;
}

interface SphereData {
  position: string | number[][] | number[];
  color: string | number[][] | number[];
  radius: string | number[] | number;
  material: string;
  resolution: number;
  scale: number;
  opacity: number;
  selecting: InteractionSettings;
  hovering: InteractionSettings;
}

const positionVec = new THREE.Vector3();
const scaleVec = new THREE.Vector3();
const matrix = new THREE.Matrix4();
const tempColor = new THREE.Color();

const HOVER_SCALE = 1.25;
const SELECTION_SCALE = 1.01;

export default function Sphere({ data, geometryKey }: { data: SphereData; geometryKey: string }) {
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
  } = data;

  const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const hoverMeshRef = useRef<THREE.Mesh | null>(null);
  const [hoveredInstanceId, setHoveredInstanceId] = useState<number | null>(null);

  const { currentFrame, roomId, clientId, selection, updateSelection, setDrawingPointerPosition, isDrawing, setDrawingIsValid } = useAppStore();
  const lastGoodFrameData = useRef<any>(null);

  const selectionSet = useMemo(() => new Set(selection || []), [selection]);
  const selectedIndices = useMemo(() => Array.from(selectionSet), [selectionSet]);

  const particleResolution = resolution || 8;
  const particleScale = scale || 1.0;

  // Determine what to fetch (skip hex colors)
  const keysToFetch = useMemo(() => {
    const keys: string[] = [];
    if (typeof positionProp === "string") keys.push(positionProp);
    if (typeof colorProp === "string" && shouldFetchAsFrameData(colorProp)) {
      keys.push(colorProp);
    }
    if (typeof radiusProp === "string") keys.push(radiusProp);
    return keys;
  }, [positionProp, colorProp, radiusProp]);

  // Fetch frame data
  const queries = useMemo(() => {
    if (!roomId) return [];
    return keysToFetch.map((key) => getFrameDataOptions(roomId, currentFrame, key));
  }, [currentFrame, roomId, keysToFetch]);

  const queryResults = useQueries({ queries });

  // Process fetched + local data into consistent format
  const { processedData, isFetching } = useMemo(() => {
    const isQueryFetching = queryResults.some((r) => r.isFetching || r.isPlaceholderData);
    const allSucceeded = queryResults.every((r) => r.isSuccess);

    if (keysToFetch.length > 0 && !allSucceeded) {
      return { isFetching: isQueryFetching, processedData: null };
    }

    const fetchedMap = new Map(queryResults.map((r, i) => [keysToFetch[i], r.data?.data]));

    let count = 0;
    let finalPositions: number[] = [];
    if (typeof positionProp === "string") {
      const remoteData = fetchedMap.get(positionProp) || [];
      count = remoteData.length / 3;
      finalPositions = remoteData;
    } else if (positionProp.length > 0 && Array.isArray(positionProp[0])) {
      count = positionProp.length;
      finalPositions = (positionProp as number[][]).flat();
    } else if (positionProp.length > 0) {
      count = 1;
      finalPositions = positionProp as number[];
    }

    if (count === 0) return { isFetching: isQueryFetching, processedData: null };

    let finalColors: number[];
    if (typeof colorProp === "string" && shouldFetchAsFrameData(colorProp)) {
      const fetchedColors = fetchedMap.get(colorProp);
      finalColors = fetchedColors ? Array.from(fetchedColors) : [];
    } else if (typeof colorProp === "string") {
      // Hex color string - convert to RGB and replicate for all particles
      const rgb = hexToRgb(colorProp);
      if (rgb) {
        finalColors = new Array(count * 3);
        for (let i = 0; i < count; i++) {
          finalColors[i * 3] = rgb[0];
          finalColors[i * 3 + 1] = rgb[1];
          finalColors[i * 3 + 2] = rgb[2];
        }
      } else {
        // Fallback to white if hex parsing fails
        throw new Error(`Invalid hex color: ${colorProp}`);
      }
    } else if (colorProp.length > 0 && Array.isArray(colorProp[0])) {
      finalColors = (colorProp as number[][]).flat();
    } else {
      finalColors = new Array(count).fill(colorProp).flat();
    }

    let finalRadii: number[];
    if (typeof radiusProp === "string") {
      finalRadii = fetchedMap.get(radiusProp) || [];
    } else if (Array.isArray(radiusProp)) {
      finalRadii = radiusProp as number[];
    } else {
      finalRadii = new Array(count).fill(radiusProp);
    }

    if (
      finalPositions.length / 3 !== count ||
      finalColors.length / 3 !== count ||
      finalRadii.length !== count
    ) {
      console.error("Sphere data arrays have inconsistent lengths.");
      return { isFetching: isQueryFetching, processedData: null };
    }

    return {
      isFetching: isQueryFetching,
      processedData: { positions: finalPositions, colors: finalColors, radii: finalRadii, count },
    };
  }, [queryResults, keysToFetch, positionProp, colorProp, radiusProp]);

  const dataToRender = processedData || lastGoodFrameData.current;
  if (processedData) lastGoodFrameData.current = processedData;

  // ⚡ Update geometry only when data changes
  useEffect(() => {
    if (!mainMeshRef.current || !dataToRender) return;

    const mainMesh = mainMeshRef.current;
    const { positions, colors, radii, count } = dataToRender;
    console.log("Rendering", count, "particles.");
    for (let i = 0; i < count; i++) {
      const i3 = i * 3;
      positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
      const r = radii[i] * particleScale;
      scaleVec.set(r, r, r);
      matrix.identity().setPosition(positionVec).scale(scaleVec);
      mainMesh.setMatrixAt(i, matrix);
      tempColor.setRGB(colors[i3], colors[i3 + 1], colors[i3 + 2]);
      mainMesh.setColorAt(i, tempColor);
    }

    mainMesh.count = count;
    mainMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    mainMesh.instanceMatrix.needsUpdate = true;
    if (mainMesh.instanceColor) mainMesh.instanceColor.needsUpdate = true;
  }, [dataToRender, particleScale]);

  useEffect(() => {
    if (!selectionMeshRef.current || !dataToRender || !selecting.enabled) return;

    const selectionMesh = selectionMeshRef.current;
    const { positions, radii, count } = dataToRender;

    selectedIndices.forEach((id, index) => {
      if (id >= count) return;
      const i3 = id * 3;
      positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
      const r = radii[id] * particleScale * SELECTION_SCALE;
      scaleVec.set(r, r, r);
      matrix.identity().setPosition(positionVec).scale(scaleVec);
      selectionMesh.setMatrixAt(index, matrix);
    });

    selectionMesh.count = selectedIndices.length;
    selectionMesh.instanceMatrix.needsUpdate = true;

  }, [selectedIndices, dataToRender, particleScale, selecting]);

  useEffect(() => {
    if (!dataToRender || isFetching) return;

    const hoverMesh = hoverMeshRef.current;
    const { positions, radii, count } = dataToRender;
    if (hoverMesh) {
      if (hovering && hoveredInstanceId !== null && hoveredInstanceId < count) {
        hoverMesh.visible = true;
        const i3 = hoveredInstanceId * 3;
        positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
        const r = radii[hoveredInstanceId] * particleScale * HOVER_SCALE;
        hoverMesh.position.copy(positionVec);
        hoverMesh.scale.set(r, r, r);
      } else {
        hoverMesh.visible = false;
      }
    }
  }, [hoveredInstanceId, dataToRender, particleScale, hovering]);


  const mainGeometry = useMemo(() => {
    return new THREE.SphereGeometry(1, particleResolution, particleResolution);
  }, [particleResolution]);


  const onClickHandler = useCallback((event: any) => {
    if (event.detail !== 1 || event.instanceId === undefined) return;
    event.stopPropagation();
    updateSelection(event.instanceId, event.shiftKey);
  }, [updateSelection]);

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
    setHoveredInstanceId(event.instanceId);
    setDrawingIsValid(true);
  }, [setDrawingIsValid]);

  const onPointerOutHandler = useCallback(() => { setHoveredInstanceId(null); setDrawingIsValid(false); }, [setDrawingIsValid]);

  if (!clientId || !roomId || !dataToRender) return null;

  return (
    <group>
      {/* Main mesh */}
      <instancedMesh
        key={dataToRender.count}
        ref={mainMeshRef}
        args={[undefined, undefined, dataToRender.count]}
        onClick={selecting.enabled ? onClickHandler : undefined}
        onPointerEnter={hovering ? onPointerEnterHandler : undefined}
        onPointerMove={hovering ? onPointerMoveHandler : undefined}
        onPointerOut={hovering ? onPointerOutHandler : undefined}
      >
        <primitive object={mainGeometry} attach="geometry" />
        {renderMaterial(material, data.opacity)}
      </instancedMesh>

      {/* Selection mesh */}
      {selecting.enabled && (
        <instancedMesh
          key={`selection-${selectedIndices.length}`}
          ref={selectionMeshRef}
          args={[undefined, undefined, selectedIndices.length]}
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

      {/* Hover mesh */}
      {hovering && (
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
    </group>
  );
}