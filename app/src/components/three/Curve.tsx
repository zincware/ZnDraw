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
  opacity: number;
}

interface CurveData {
  position: string | number[][];
  color: string;
  material: "LineBasicMaterial" | "LineDashedMaterial";
  divisions: number;
  variant: "CatmullRomCurve3";
  thickness: number;
  marker: MarkerData | null;
  virtual_marker: MarkerData | null;
}

/**
 * A component to render a 3D curve using @react-three/drei,
 * with optional markers that can be interactively moved.
 */
export default function Curve({ data, geometryKey }: { data: CurveData; geometryKey: string }) {
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
  const virtualMarkerRefs = useRef<(THREE.Mesh | null)[]>([]);

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
    virtualMarkerRefs.current = virtualMarkerRefs.current.slice(0, Math.max(0, sourceCurvePoints.length - 1));
  }, [sourceCurvePoints]);

  // --- Calculate virtual marker positions ON the curve between each pair of control points ---
  const virtualMarkerPositions = useMemo(() => {
    if (!interactivePoints || interactivePoints.length < 2) return [];
    const curve = new THREE.CatmullRomCurve3(interactivePoints);
    const positions: THREE.Vector3[] = [];
    
    // Get all points on the curve for analysis
    const curvePoints = curve.getPoints(data.divisions * interactivePoints.length);
    
    // For each pair of consecutive control points, find the midpoint on the curve
    for (let i = 0; i < interactivePoints.length - 1; i++) {
      const controlPoint1 = interactivePoints[i];
      const controlPoint2 = interactivePoints[i + 1];
      
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
      positions.push(curvePoints[midIndex].clone());
    }
    return positions;
  }, [interactivePoints, data.divisions]);

  // --- Handler to convert virtual marker to actual marker ---
  const handleVirtualMarkerClick = (insertIndex: number) => {
    if (!interactivePoints) return;
    const newPoints = [...interactivePoints];
    newPoints.splice(insertIndex + 1, 0, virtualMarkerPositions[insertIndex].clone());
    setInteractivePoints(newPoints);
    // Set the newly created point as selected
    setSelectedIndex(insertIndex + 1);
  };

  // --- Handler to send updated geometry data to backend ---
  const handleMouseUp = async () => {
    if (!roomId || !interactivePoints) return;

    // Convert Vector3 array to array of tuples
    const positions = interactivePoints.map(point => [point.x, point.y, point.z]);

    // Build the complete data object matching the backend format
    const geometryData = {
      position: positions,
      color: data.color,
      material: data.material,
      variant: data.variant,
      divisions: data.divisions,
      thickness: data.thickness,
      marker: data.marker,
      virtual_marker: data.virtual_marker,
    };

    try {
      const response = await fetch(`/api/rooms/${roomId}/geometries`, {
        method: "PUT",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          key: geometryKey,
          data: geometryData,
          type: "Curve",
        }),
      });

      if (!response.ok) {
        console.error("Failed to update geometry:", await response.text());
      }
    } catch (error) {
      console.error("Error updating geometry:", error);
    }
  };

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
            key={`marker-${index}`}
            ref={(el) => {
              markerRefs.current[index] = el;
            }}
            position={point}
            args={[marker.size]}
            onClick={(e) => {
              e.stopPropagation();
              setSelectedIndex(index);
            }}
          >
            <meshBasicMaterial 
              color={marker.color || color} 
              opacity={marker.opacity}
              transparent={marker.opacity < 1}
            />
          </Dodecahedron>
        ))}

      {/* Virtual Markers (midpoints between consecutive markers) */}
      {data.virtual_marker && virtualMarkerPositions.map((position, index) => {
        const virtualMarker = data.virtual_marker!;
        return (
          <Dodecahedron
            key={`virtual-${index}`}
            ref={(el) => {
              virtualMarkerRefs.current[index] = el;
            }}
            position={position}
            args={[virtualMarker.size]}
            onClick={(e) => {
              e.stopPropagation();
              handleVirtualMarkerClick(index);
            }}
          >
            <meshBasicMaterial 
              color={virtualMarker.color || color} 
              opacity={virtualMarker.opacity}
              transparent={virtualMarker.opacity < 1}
            />
          </Dodecahedron>
        );
      })}

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
          onMouseUp={handleMouseUp}
          onPointerMissed={() => setSelectedIndex(null)}
        />
      )}
    </group>
  );
}