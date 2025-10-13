import { useMemo, useRef, useState, useEffect, useCallback } from "react";
import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { useTheme } from "@mui/material/styles";
import {
  Dodecahedron,
  TransformControls,
} from "@react-three/drei";
import { Line } from "@react-three/drei";
import { useAppStore } from "../../store";
import { getFrames, createGeometry } from "../../myapi/client";
import { debounce } from "lodash";

interface MarkerData {
  enabled: boolean;
  size: number;
  color: string;
  opacity: number;
}

interface CurveData {
  position: string | number[][];
  color: string;
  material: "LineBasicMaterial" | "LineDashedMaterial";
  divisions: number;
  variant: "CatmullRomCurve3";
  thickness: number;
  marker: MarkerData;
  virtual_marker: MarkerData;
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

  const theme = useTheme();

  const { currentFrame, roomId, isDrawing, drawingIsValid, drawingPointerPosition, setIsDrawing, clientId, setGeometryFetching, removeGeometryFetching, setCurveLength } = useAppStore();
  const [markerPositions, setMarkerPositions] = useState<THREE.Vector3[]>([]);
  const [lineSegments, setLineSegments] = useState<THREE.Vector3[]>([]);
  const [virtualMarkerPositions, setVirtualMarkerPositions] = useState<THREE.Vector3[]>([]);
  const [selectedIndex, setSelectedIndex] = useState<number | null>(null);

  const [lastUpdateSource, setLastUpdateSource] = useState<'remote' | 'local' | null>(null);

  const markerRefs = useRef<(THREE.Mesh | null)[]>([]);

  // Helper function to resolve color based on theme
  const getLineColor = useCallback((colorValue: string) => {
    return colorValue === "default"
      ? theme.palette.text.primary
      : colorValue;
  }, [theme]);

  const getMarkerColor = useCallback((colorValue: string) => {
    return colorValue === "default"
      ? theme.palette.secondary.main
      : colorValue;
  }, [theme]);

  // --- Data Fetching ---
  // Simple conditional query for position data
  const { data: positionQueryData, isFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, positionProp],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [positionProp as string], signal),
    enabled: !!roomId && !!clientId && typeof positionProp === "string",
    placeholderData: keepPreviousData,
  });

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

  useEffect(() => {
    if (typeof positionProp === "string") {
      // Fetched from server
      if (isFetching) return;

      const points = positionQueryData?.[positionProp];
      if (!points) {
        setMarkerPositions([]);
        setLastUpdateSource('remote');
        return;
      }

      const vecPoints = [];
      for (let i = 0; i < points.length; i += 3) {
        vecPoints.push(new THREE.Vector3(points[i], points[i + 1], points[i + 2]));
      }
      setMarkerPositions(vecPoints);
      setLastUpdateSource('remote');
    } else {
      // Static data
      const manualPoints = positionProp as number[][];
      if (manualPoints === undefined) {
        console.error("Manual points are undefined");
        return;
      }
      const vecPoints = manualPoints.map(p => new THREE.Vector3(p[0], p[1], p[2]));
      setMarkerPositions(vecPoints);
      setLastUpdateSource('remote');
    }
  }, [positionQueryData, isFetching, positionProp]);

  useEffect(() => {
    const allPoints: THREE.Vector3[] = [...markerPositions];
    if (!data) return;
    if (drawingPointerPosition !== null && isDrawing) {
      allPoints.push(drawingPointerPosition);
    }
    if (allPoints.length < 2) {
      setLineSegments([]);
      setVirtualMarkerPositions([]);
      setCurveLength(0);
      return;
    }
    const curve = new THREE.CatmullRomCurve3(allPoints);
    const curvePoints = curve.getPoints(data.divisions * allPoints.length);
    setLineSegments(curvePoints);

    setCurveLength(curve.getLength());

    const _virtualMarkerPositions: THREE.Vector3[] = [];
    // virtual markers
    for (let i = 0; i < allPoints.length - 1; i++) {
      const controlPoint1 = allPoints[i];
      const controlPoint2 = allPoints[i + 1];

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
      _virtualMarkerPositions.push(curvePoints[midIndex].clone());
    }
    setVirtualMarkerPositions(_virtualMarkerPositions);

  }, [markerPositions, drawingPointerPosition, isDrawing, data.divisions, setCurveLength]);

  const handleVirtualMarkerClick = useCallback((insertIndex: number) => {
    const newPoints = [...markerPositions];
    newPoints.splice(insertIndex + 1, 0, virtualMarkerPositions[insertIndex].clone());
    setMarkerPositions(newPoints);
    // We need to select *after* the new marker has been added
    // Delay selection to the next tick so refs update
    setTimeout(() => {
      setSelectedIndex(insertIndex + 1);
    }, 0);
  }, [markerPositions, virtualMarkerPositions]);

  const handleDrawingMarkerClick = useCallback(() => {
    if (drawingIsValid && drawingPointerPosition) {
      const newPoints = [...markerPositions];
      newPoints.push(drawingPointerPosition.clone());
      setLastUpdateSource('local');
      setMarkerPositions(newPoints);
    } else {
      setIsDrawing(false);
    }
  }, [drawingIsValid, drawingPointerPosition, markerPositions, setIsDrawing]);

  const selectedMesh = useMemo(() => {
    return selectedIndex !== null ? markerRefs.current[selectedIndex] : null;
  }, [selectedIndex]);

  // attach event handler to backspace/delete key to remove selected marker
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === "Backspace" || e.key === "Delete") {
        // Use functional state update to avoid stale closure
        setSelectedIndex(null);

        // We need to remove the transform controls first,
        // because the on-change handler would recreate the marker
        setTimeout(() => {
          setLastUpdateSource('local');
          setMarkerPositions(prev => {
            if (selectedIndex === null) return prev;
            const newPoints = [...prev];
            newPoints.splice(selectedIndex, 1);
            return newPoints;
          });
        }, 0);
      }
    };

    window.addEventListener("keydown", handleKeyDown);
    return () => {
      window.removeEventListener("keydown", handleKeyDown);
    };
    // 👇 only run once; don't depend on selectedIndex or markerPositions
  }, [selectedIndex]);

  const persistMarkerPositions = useCallback(async () => {
    if (!roomId || !markerPositions) return;

    // Convert Vector3 array to array of tuples
    const positions = markerPositions.map(point => [point.x, point.y, point.z]);

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
      await createGeometry(roomId, clientId, geometryKey, "Curve", geometryData);
    } catch (error) {
      console.error("Error updating geometry:", error);
    }
  }, [roomId, markerPositions, geometryKey, data]);

  useEffect(() => {
    if (lastUpdateSource !== 'local') return;
    
    const debouncedPersist = debounce(persistMarkerPositions, 500);
    debouncedPersist();
    return () => {
      debouncedPersist.cancel();
    };
  }, [markerPositions, persistMarkerPositions, lastUpdateSource]);

  const onTransformChange = useCallback(() => {
    if (selectedIndex !== null && markerRefs.current[selectedIndex]) {
      setMarkerPositions(prev => {
        const newPoints = [...prev];
        newPoints[selectedIndex] = markerRefs.current[selectedIndex]!.position.clone();
        return newPoints;
      });
      setLastUpdateSource('local');
    }
  }, [selectedIndex]);

  if (!roomId) return null;
  if (!marker) return null;

  // --- Render ---
  return (
    <group>
      {/* Transform Controls */}
      {selectedMesh && selectedIndex !== null &&
        <TransformControls
          object={selectedMesh}
          onPointerMissed={() => setSelectedIndex(null)}
          mode="translate"
          onChange={onTransformChange}
        />
      }

      {/* Render curve line */}
      {lineSegments.length >= 2 && (
        <Line
          points={lineSegments}
          color={(isDrawing && !drawingIsValid) ? "#FF0000" : getLineColor(color)}
          lineWidth={thickness}
          dashed={material === "LineDashedMaterial"}
        />
      )}

      {/* Markers */}
      {marker.enabled && markerPositions &&
        markerPositions.map((point, index) => (
          <Dodecahedron
            key={`marker-${index}`}
            position={point}
            args={[marker.size]}
            ref={(el) => {
              markerRefs.current[index] = el;
            }}
            onClick={(e) => {
              e.stopPropagation();
              setSelectedIndex(index);
            }}
          >
            <meshBasicMaterial
              color={(isDrawing && !drawingIsValid) ? "#FF0000" : getMarkerColor(marker.color)}
              opacity={marker.opacity}
              transparent={marker.opacity < 1}
            />
          </Dodecahedron>
        ))}

      {/* Virtual Markers */}
      {data.virtual_marker.enabled && virtualMarkerPositions.map((position, index) => {
        const virtualMarker = data.virtual_marker!;
        return (
          <Dodecahedron
            key={`virtual-${index}`}
            position={position}
            args={[virtualMarker.size]}
            onClick={() => handleVirtualMarkerClick(index)}
          >
            <meshBasicMaterial
              color={getMarkerColor(virtualMarker.color)}
              opacity={virtualMarker.opacity}
              transparent={virtualMarker.opacity < 1}
            />
          </Dodecahedron>
        );
      })}

      {/* Drawing Marker */}
      {isDrawing && drawingPointerPosition && marker.enabled && (
        <Dodecahedron
          key={`virtual-drawing`}
          position={drawingPointerPosition}
          args={[marker.size]}
          onClick={handleDrawingMarkerClick}
        >
          <meshBasicMaterial
            color={(isDrawing && !drawingIsValid) ? "#FF0000" : getMarkerColor(marker.color)}
            opacity={marker.opacity}
            transparent={marker.opacity < 1}
          />
        </Dodecahedron>
      )}
    </group>
  );
}