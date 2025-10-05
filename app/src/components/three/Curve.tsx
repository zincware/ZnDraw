import { useMemo, useRef, useState, useEffect } from "react";
import * as THREE from "three";
import { useQueries } from "@tanstack/react-query";
import {
  Dodecahedron,
  TransformControls,
} from "@react-three/drei";
import { Line } from "@react-three/drei";
import { getFrameDataOptions } from "../../hooks/useTrajectoryData";
import { useAppStore } from "../../store";

interface MarkerData {
  size: number;
  color: string | null;
}

interface CurveData {
  position: string | number[][];
  color: string;
  material: "LineBasicMaterial" | "LineDashedMaterial";
  divisions: number;
  variant: "CatmullRomCurve3";
  thickness: number;
  marker: MarkerData | null;
}

/**
 * A component to render a 3D curve using @react-three/drei,
 * with optional markers that can be interactively moved.
 */
export default function Curve({ data }: { data: CurveData }) {
  const {
    position: positionProp,
    color,
    material,
    thickness,
    marker,
  } = data;

  const { currentFrame, roomId } = useAppStore();
  const lastGoodFrameData = useRef<{ points: number[] } | null>(null);

  // --- State for Interactivity ---
  const [interactivePoints, setInteractivePoints] = useState<THREE.Vector3[] | null>(null);
  const [selectedIndex, setSelectedIndex] = useState<number | null>(null);
  const markerRefs = useRef<(THREE.Mesh | null)[]>([]);

  // --- Data Fetching ---
  const shouldFetchPosition = typeof positionProp === "string";

  const queries = useMemo(() => {
    if (!roomId || !shouldFetchPosition) return [];
    return [getFrameDataOptions(roomId, currentFrame, positionProp as string)];
  }, [currentFrame, roomId, positionProp, shouldFetchPosition]);

  const queryResults = useQueries({ queries });

  const processedData = useMemo(() => {
    if (shouldFetchPosition) {
      const result = queryResults[0];
      if (!result || !result.isSuccess) return null;
      const points = result.data?.data || [];
      return points.length < 6 ? null : { points };
    } else {
      const manualPoints = positionProp as number[][];
      if (manualPoints.length < 2) return null;
      const points = manualPoints.flat();
      return { points };
    }
  }, [queryResults, positionProp, shouldFetchPosition]);

  // Keep last good frame data so the line doesn’t flicker between frames
  const dataToRender = processedData || lastGoodFrameData.current;
  if (processedData) lastGoodFrameData.current = processedData;

  // --- Stable memoization of 3D points ---
  // Converts flat array → Vector3 array, but only when the underlying numeric values actually change.
  const sourceCurvePoints = useMemo(() => {
    if (!dataToRender?.points) return null;

    // Create a unique key based on numeric content.
    const key = dataToRender.points.join(",");
    const vectors: THREE.Vector3[] = [];

    for (let i = 0; i < dataToRender.points.length; i += 3) {
      vectors.push(
        new THREE.Vector3(
          dataToRender.points[i],
          dataToRender.points[i + 1],
          dataToRender.points[i + 2]
        )
      );
    }
    return vectors.length < 2 ? null : vectors;
  }, [dataToRender?.points?.join(",")]);

  // --- Sync effect ---
  // Updates interactivity state only when actual data changes (no infinite loops!)
  useEffect(() => {
    if (!sourceCurvePoints) return;
    setInteractivePoints(sourceCurvePoints);
    markerRefs.current = markerRefs.current.slice(0, sourceCurvePoints.length);
  }, [sourceCurvePoints]);

  // --- Early exit if not ready ---
  if (!roomId || !interactivePoints) return null;

  const selectedMesh = selectedIndex !== null ? markerRefs.current[selectedIndex] : null;

  // --- Generate curve points ---
  const curve = new THREE.CatmullRomCurve3(interactivePoints);
  const linePoints = curve.getPoints(data.divisions * interactivePoints.length);

  // --- Render ---
  return (
    <group>
      <Line
        points={linePoints}
        color={color}
        lineWidth={thickness}
        dashed={material === "LineDashedMaterial"}
      />

      {/* Optional Markers */}
      {marker &&
        interactivePoints.map((point, index) => (
          <Dodecahedron
            key={index}
            ref={(el) => (markerRefs.current[index] = el)}
            position={point}
            args={[marker.size]}
            onClick={(e) => {
              e.stopPropagation();
              setSelectedIndex(index);
            }}
          >
            <meshBasicMaterial color={marker.color || color} />
          </Dodecahedron>
        ))}

      {/* Transform Controls for interactivity */}
      {selectedMesh && (
        <TransformControls
          object={selectedMesh}
          onChange={() => {
            if (selectedIndex !== null && markerRefs.current[selectedIndex]) {
              const newPoints = [...interactivePoints];
              newPoints[selectedIndex] =
                markerRefs.current[selectedIndex]!.position.clone();
              setInteractivePoints(newPoints);
            }
          }}
          onPointerMissed={() => setSelectedIndex(null)}
        />
      )}
    </group>
  );
}