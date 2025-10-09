import * as THREE from "three";
import { useAppStore } from "../../store";
import { useRef, useMemo, useEffect, useState } from "react";
import { BufferGeometryUtils } from "three/examples/jsm/Addons.js";
import { renderMaterial } from "./materials";
import { getFrames } from "../../myapi/client";
import { decode } from "@msgpack/msgpack";
import { useQuery } from '@tanstack/react-query';
import { shouldFetchAsFrameData, hexToRgb } from "../../utils/colorUtils";

const numpyDtypeToTypedArray = {
  float32: Float32Array,
  float64: Float64Array,
  int8: Int8Array,
  int16: Int16Array,
  int32: Int32Array,
  uint8: Uint8Array,
  uint16: Uint16Array,
  uint32: Uint32Array,
};


// Props interface for the dynamic keys
type StaticValue = number | number[] | number[][];
type DataProp = string | StaticValue;

interface ArrowProps {
  position: DataProp;
  direction: DataProp;
  color: DataProp;
  radius: DataProp;
  scale: DataProp;
  material: string;
  geometryKey: string;
}

// Reusable THREE objects
const _startVec = new THREE.Vector3();
// No need for _endVec if we use direction directly
const _direction = new THREE.Vector3();
// The default orientation of the arrow geometry is along the +Y axis
const _arrowUp = new THREE.Vector3(0, 1, 0);
const _quaternion = new THREE.Quaternion();
const _scaleVec = new THREE.Vector3();
const _matrix = new THREE.Matrix4();
const _color = new THREE.Color();

function createArrowMesh() {
  const cylinderRadius = 0.04;
  const cylinderHeight = 0.6;
  const coneRadius = 0.1;
  const coneHeight = 0.4;

  const cylinderGeometry = new THREE.CylinderGeometry(
    cylinderRadius,
    cylinderRadius,
    cylinderHeight,
    32,
  );
  const coneGeometry = new THREE.ConeGeometry(coneRadius, coneHeight, 32);

  // Position geometries so the base is at (0,0,0) and it extends up to a total height of 1.0
  cylinderGeometry.translate(0, cylinderHeight / 2, 0);
  coneGeometry.translate(0, cylinderHeight + coneHeight / 2, 0);

  const arrowGeometry = BufferGeometryUtils.mergeGeometries([
    cylinderGeometry,
    coneGeometry,
  ]);

  // This geometry now has a total height of 1.0 and points along the +Y axis.
  return arrowGeometry;
}

export default function Arrow({ data, geometryKey }: { data: ArrowProps; geometryKey: string }) {
  const { position, direction, color, radius, scale, material } = data;
  const instancedMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const { currentFrame, roomId, clientId, selection } = useAppStore();
  const [fetchedData, setFetchedData] = useState<Record<string, any>>({});

  // Determine which keys need fetching (exclude hex colors)Cwi
  const keysToFetch = useMemo(() => {
    const keys: Record<string, string> = {};
    if (typeof position === "string") keys.position = position;
    if (typeof direction === "string") keys.direction = direction;
    if (typeof color === "string" && shouldFetchAsFrameData(color)) keys.color = color;
    if (typeof radius === "string") keys.radius = radius;
    if (typeof scale === "string") keys.scale = scale;
    console.log("Arrow keysToFetch:", keys);
    return keys;
  }, [position, direction, color, radius, scale]);

  // Create queries for each dynamic property
  const { data: positionData, isFetching: isPositionFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, keysToFetch.position],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId, currentFrame, [keysToFetch.position], signal),
    enabled: !!keysToFetch.position,
  });

  const { data: directionData, isFetching: isDirectionFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, keysToFetch.direction],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId, currentFrame, [keysToFetch.direction], signal),
    enabled: !!keysToFetch.direction,
  });

  const { data: colorData, isFetching: isColorFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, keysToFetch.color],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId, currentFrame, [keysToFetch.color], signal),
    enabled: !!keysToFetch.color,
  });

  const { data: radiusData, isFetching: isRadiusFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, keysToFetch.radius],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId, currentFrame, [keysToFetch.radius], signal),
    enabled: !!keysToFetch.radius,
  });

  const { data: scaleData, isFetching: isScaleFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, keysToFetch.scale],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId, currentFrame, [keysToFetch.scale], signal),
    enabled: !!keysToFetch.scale,
  });

  // Check if any enabled query is still fetching
  const isFetching =
    (keysToFetch.position && isPositionFetching) ||
    (keysToFetch.direction && isDirectionFetching) ||
    (keysToFetch.color && isColorFetching) ||
    (keysToFetch.radius && isRadiusFetching) ||
    (keysToFetch.scale && isScaleFetching);

  // Clean up fetchedData when keysToFetch changes (e.g., switching from data key to static value)
  useEffect(() => {
    setFetchedData((prev) => {
      const cleanedData: Record<string, any> = {};
      // Only keep data that's still in keysToFetch
      for (const key of Object.keys(keysToFetch)) {
        if (prev[key]) {
          cleanedData[key] = prev[key];
        }
      }
      return cleanedData;
    });
  }, [keysToFetch]);

  // Decode position data
  useEffect(() => {
    if (positionData && keysToFetch.position) {
      const decodedData = decode(positionData)[0][keysToFetch.position] as any;
      const TypedArray = numpyDtypeToTypedArray[decodedData.dtype as keyof typeof numpyDtypeToTypedArray];
      const dataArray = new TypedArray(decodedData.data.slice().buffer);

      setFetchedData((prev) => ({
        ...prev,
        position: dataArray,
      }));
    }
  }, [positionData, keysToFetch.position]);

  // Decode direction data
  useEffect(() => {
    if (directionData && keysToFetch.direction) {
      const decodedData = decode(directionData)[0][keysToFetch.direction] as any;
      const TypedArray = numpyDtypeToTypedArray[decodedData.dtype as keyof typeof numpyDtypeToTypedArray];
      const dataArray = new TypedArray(decodedData.data.slice().buffer);

      setFetchedData((prev) => ({
        ...prev,
        direction: dataArray,
      }));
    }
  }, [directionData, keysToFetch.direction]);

  // Decode color data
  useEffect(() => {
    if (colorData && keysToFetch.color) {
      const decodedData = decode(colorData)[0][keysToFetch.color] as any;
      const TypedArray = numpyDtypeToTypedArray[decodedData.dtype as keyof typeof numpyDtypeToTypedArray];
      const dataArray = new TypedArray(decodedData.data.slice().buffer);

      setFetchedData((prev) => ({
        ...prev,
        color: dataArray,
      }));
    }
  }, [colorData, keysToFetch.color]);

  // Decode radius data
  useEffect(() => {
    if (radiusData && keysToFetch.radius) {
      const decodedData = decode(radiusData)[0][keysToFetch.radius] as any;
      const TypedArray = numpyDtypeToTypedArray[decodedData.dtype as keyof typeof numpyDtypeToTypedArray];
      const dataArray = new TypedArray(decodedData.data.slice().buffer);

      setFetchedData((prev) => ({
        ...prev,
        radius: dataArray,
      }));
    }
  }, [radiusData, keysToFetch.radius]);

  // Decode scale data
  useEffect(() => {
    if (scaleData && keysToFetch.scale) {
      const decodedData = decode(scaleData)[0][keysToFetch.scale] as any;
      const TypedArray = numpyDtypeToTypedArray[decodedData.dtype as keyof typeof numpyDtypeToTypedArray];
      const dataArray = new TypedArray(decodedData.data.slice().buffer);

      setFetchedData((prev) => ({
        ...prev,
        scale: dataArray,
      }));
    }
  }, [scaleData, keysToFetch.scale]);

  // Process fetched + static data into final format
  const processedData = useMemo(() => {
    console.log("Processing Arrow data with fetchedData:", fetchedData);
    try {
      // Determine count from fetched or static data
      let finalCount = 0;
      let finalPositions: number[] = [];
      let finalDirections: number[] = [];
      let finalColors: number[] = [];
      let finalRadii: number[] = [];
      let finalScales: number[] = [];

      // Process position
      if (fetchedData.position) {
        finalCount = fetchedData.position.length / 3;
        finalPositions = Array.from(fetchedData.position);
      } else if (typeof position !== "string") {
        if (Array.isArray(position[0])) {
          finalCount = position.length;
          finalPositions = (position as number[][]).flat();
        } else {
          finalCount = 1;
          finalPositions = position as number[];
        }
      }

      if (finalCount === 0) return null;

      // Process direction
      if (fetchedData.direction) {
        finalDirections = Array.from(fetchedData.direction);
      } else if (typeof direction !== "string") {
        if (Array.isArray(direction[0])) {
          finalDirections = (direction as number[][]).flat();
        } else if (Array.isArray(direction)) {
          finalDirections = new Array(finalCount * 3);
          for (let i = 0; i < finalCount; i++) {
            finalDirections[i * 3] = (direction as number[])[0];
            finalDirections[i * 3 + 1] = (direction as number[])[1];
            finalDirections[i * 3 + 2] = (direction as number[])[2];
          }
        }
      }

      // Process color
      if (fetchedData.color) {
        finalColors = Array.from(fetchedData.color);
      } else if (typeof color === "string") {
        // Hex color
        const rgb = hexToRgb(color);
        if (rgb) {
          finalColors = new Array(finalCount * 3);
          for (let i = 0; i < finalCount; i++) {
            finalColors[i * 3] = rgb[0];
            finalColors[i * 3 + 1] = rgb[1];
            finalColors[i * 3 + 2] = rgb[2];
          }
        } else {
          throw new Error(`Invalid hex color: ${color}`);
        }
      } else if (Array.isArray(color[0])) {
        finalColors = (color as number[][]).flat();
      } else if (Array.isArray(color)) {
        finalColors = new Array(finalCount * 3);
        for (let i = 0; i < finalCount; i++) {
          finalColors[i * 3] = (color as number[])[0];
          finalColors[i * 3 + 1] = (color as number[])[1];
          finalColors[i * 3 + 2] = (color as number[])[2];
        }
      }

      // Process radius
      if (fetchedData.radius) {
        finalRadii = Array.from(fetchedData.radius);
      } else if (typeof radius !== "string") {
        if (Array.isArray(radius)) {
          finalRadii = radius as number[];
        } else {
          finalRadii = new Array(finalCount).fill(radius);
        }
      }

      // Process scale
      if (fetchedData.scale) {
        finalScales = Array.from(fetchedData.scale);
      } else if (typeof scale !== "string") {
        if (Array.isArray(scale)) {
          finalScales = scale as number[];
        } else {
          finalScales = new Array(finalCount).fill(scale);
        }
      }

      // Validate lengths
      if (
        finalPositions.length / 3 !== finalCount ||
        finalDirections.length / 3 !== finalCount ||
        finalColors.length / 3 !== finalCount ||
        finalRadii.length !== finalCount ||
        finalScales.length !== finalCount
      ) {
        console.error("Arrow data arrays have inconsistent lengths");
        return null;
      }

      return {
        positions: finalPositions,
        directions: finalDirections,
        colors: finalColors,
        radii: finalRadii,
        scales: finalScales,
        count: finalCount,
      };
    } catch (error) {
      console.error("Error processing Arrow data:", error);
      return null;
    }
  }, [fetchedData, position, direction, color, radius, scale]);


  // Update arrow instances when data changes
  useEffect(() => {
    if (!processedData || !instancedMeshRef.current) return;

    const mesh = instancedMeshRef.current;
    const { positions, directions, colors, radii, scales, count: dataCount } = processedData;

    console.log("Updating", dataCount, "arrows");

    const selectionSet = selection ? new Set(selection) : null;

    for (let i = 0; i < dataCount; i++) {
      const i3 = i * 3;

      // Set the starting position of the arrow
      _startVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);

      // Set the direction vector
      _direction.set(directions[i3], directions[i3 + 1], directions[i3 + 2]);

      // Calculate the length of the arrow from the direction vector's magnitude
      const arrowLength = _direction.length() * scales[i];

      // Set the scale (Y-axis is height/length, X and Z are width/radius)
      const arrowRadius = radii[i];
      _scaleVec.set(arrowRadius, arrowLength, arrowRadius);

      // Calculate the rotation quaternion
      if (arrowLength > 0.0001) {
        _quaternion.setFromUnitVectors(_arrowUp, _direction.normalize());
      } else {
        _quaternion.identity();
      }

      // Compose the final transformation matrix
      _matrix.compose(_startVec, _quaternion, _scaleVec);
      mesh.setMatrixAt(i, _matrix);

      // Set color (with selection highlight)
      if (selectionSet && selectionSet.has(i)) {
        _color.setRGB(1.0, 0.75, 0.8); // Pink for selection
      } else {
        _color.setRGB(colors[i3], colors[i3 + 1], colors[i3 + 2]);
      }
      mesh.setColorAt(i, _color);
    }

    mesh.instanceMatrix.needsUpdate = true;
    if (mesh.instanceColor) {
      mesh.instanceColor.needsUpdate = true;
    }
  }, [processedData, selection]);

  const geometry = useMemo(createArrowMesh, []);

  // Return null if no data is available yet
  if (!processedData || !clientId || !roomId) {
    return null;
  }

  return (
    <instancedMesh
      key={processedData.count}
      ref={instancedMeshRef}
      args={[geometry, undefined, processedData.count]}
    >
      {renderMaterial(material)}
    </instancedMesh>
  );
}
