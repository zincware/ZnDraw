import { useMemo, useRef } from "react";
import * as THREE from "three"; // Still needed for THREE.Vector3
import { useQueries } from "@tanstack/react-query";
import { CatmullRomLine, Dodecahedron } from "@react-three/drei";
import { getFrameDataOptions } from "../../hooks/useTrajectoryData"; // Assuming this hook exists
import { useAppStore } from "../../store"; // Assuming this store exists

interface MarkerData {
    size: number;
    color: string | null;
}

// Define the structure for the curve's properties
interface CurveData {
  position: string | number[][];
  color: string;
  material: "LineBasicMaterial" | "LineDashedMaterial";
  divisions: number; // Note: CatmullRomLine handles divisions internally. This prop is no longer directly used for geometry.
  variant: "CatmullRomCurve3";
  thickness: number;
  marker: MarkerData | null;
}

/**
 * A component to render a 3D curve using @react-three/drei,
 * with optional markers at each control point.
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

  // --- Data Fetching ---

  const shouldFetchPosition = typeof positionProp === "string";

  const queries = useMemo(() => {
    if (!roomId || !shouldFetchPosition) return [];
    return [getFrameDataOptions(roomId, currentFrame, positionProp as string)];
  }, [currentFrame, roomId, positionProp, shouldFetchPosition]);

  const queryResults = useQueries({ queries });

  // --- Data Processing ---

  const { processedData } = useMemo(() => {
    const isQueryFetching = queryResults.some((r) => r.isFetching || r.isPlaceholderData);

    if (shouldFetchPosition) {
      const result = queryResults[0];
      if (!result || !result.isSuccess) {
        return { processedData: null };
      }
      const points = result.data?.data || [];
      if (points.length < 6) return { processedData: null };
      return { processedData: { points } };
    } else {
      if (positionProp.length < 2) {
         return { processedData: null };
      }
      const points = (positionProp as number[][]).flat();
      return { processedData: { points } };
    }
  }, [queryResults, positionProp, shouldFetchPosition]);

  const dataToRender = processedData || lastGoodFrameData.current;
  if (processedData) {
    lastGoodFrameData.current = processedData;
  }

  // --- Points Generation for Drei ---

  // Memoize the array of THREE.Vector3 points, which will serve as the control points.
  const curvePoints = useMemo(() => {
    if (!dataToRender) return null;

    const { points: flatPoints } = dataToRender;
    const vectors = [];
    for (let i = 0; i < flatPoints.length; i += 3) {
      vectors.push(new THREE.Vector3(flatPoints[i], flatPoints[i + 1], flatPoints[i + 2]));
    }

    // A curve needs at least two control points.
    if (vectors.length < 2) return null;
    
    return vectors;
  }, [dataToRender]);


  // --- Rendering ---

  // If there's no data to render, render nothing.
  if (!roomId || !curvePoints) return null;

  return (
    // Use a key to ensure the component remounts if the number of points changes.
    // A group is used to contain both the line and its markers.
    <group key={dataToRender?.points.length}>
      <CatmullRomLine
        points={curvePoints}
        color={color}
        lineWidth={thickness}
        dashed={material === "LineDashedMaterial"}
        // Other CatmullRomLine props like `tension` or `dashScale` could be added here
      />

      {/* Conditionally render markers at each control point */}
      {marker && curvePoints.map((point, index) => (
        <Dodecahedron key={index} position={point} args={[marker.size]}>
            <meshBasicMaterial color={marker.color || color} />
        </Dodecahedron>
      ))}
    </group>
  );
}
