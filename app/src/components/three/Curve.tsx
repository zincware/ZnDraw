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
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

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
  const { geometryDefaults } = useAppStore();

  // Merge with defaults from Pydantic (single source of truth)
  // Memoize the entire result to prevent re-creation on every render
  const fullData = useMemo(
    () => getGeometryWithDefaults<CurveData>(data, "Curve", geometryDefaults),
    [JSON.stringify(data), JSON.stringify(geometryDefaults.Curve)]
  );

  // Memoize object/array references to prevent infinite re-renders
  const positionProp = useMemo(() => fullData.position, [
    typeof fullData.position === 'string' ? fullData.position : JSON.stringify(fullData.position)
  ]);
  const marker = useMemo(() => fullData.marker, [JSON.stringify(fullData.marker)]);
  const virtual_marker = useMemo(() => fullData.virtual_marker, [JSON.stringify(fullData.virtual_marker)]);

  const {
    color,
    material,
    thickness,
    divisions,
    variant,
  } = fullData;

  const theme = useTheme();

  const { currentFrame, roomId, isDrawing, drawingIsValid, drawingPointerPosition, setIsDrawing, clientId, setGeometryFetching, removeGeometryFetching, setCurveLength, activeCurveForDrawing, registerCurveRef, unregisterCurveRef } = useAppStore();
  const [markerPositions, setMarkerPositions] = useState<THREE.Vector3[]>([]);
  const [lineSegments, setLineSegments] = useState<THREE.Vector3[]>([]);
  const [virtualMarkerPositions, setVirtualMarkerPositions] = useState<THREE.Vector3[]>([]);
  const [selectedIndex, setSelectedIndex] = useState<number | null>(null);

  const [lastUpdateSource, setLastUpdateSource] = useState<'remote' | 'local' | null>(null);

  const markerRefs = useRef<(THREE.Mesh | null)[]>([]);
  const curveRef = useRef<THREE.CatmullRomCurve3 | null>(null);

  // Check if this curve is the active drawing target
  const isActiveDrawingTarget = activeCurveForDrawing === geometryKey;

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

  // Clean up fetching state and curve ref on unmount
  useEffect(() => {
    return () => {
      removeGeometryFetching(geometryKey);
      unregisterCurveRef(geometryKey);
    };
  }, [geometryKey, removeGeometryFetching, unregisterCurveRef]);

  useEffect(() => {
    if (typeof positionProp === "string") {
      // Fetched from server
      if (isFetching) return;

      const points = positionQueryData?.[positionProp];
      if (!points) {
        setMarkerPositions(prev => prev.length === 0 ? prev : []);
        setLastUpdateSource('remote');
        return;
      }

      const vecPoints = [];
      for (let i = 0; i < points.length; i += 3) {
        vecPoints.push(new THREE.Vector3(
          Number(points[i]),
          Number(points[i + 1]),
          Number(points[i + 2])
        ));
      }
      setMarkerPositions(vecPoints);
      setLastUpdateSource('remote');
    } else {
      // Static data
      const manualPoints = positionProp as number[][];
      if (!manualPoints || manualPoints.length === 0) {
        // Don't update if already empty - prevents infinite loop
        setMarkerPositions(prev => prev.length === 0 ? prev : []);
        setLastUpdateSource('remote');
        return;
      }
      const vecPoints = manualPoints.map(p => new THREE.Vector3(p[0], p[1], p[2]));
      setMarkerPositions(vecPoints);
      setLastUpdateSource('remote');
    }
  }, [positionQueryData, isFetching, positionProp]);

  useEffect(() => {
    const allPoints: THREE.Vector3[] = [...markerPositions];
    if (!fullData) return;
    // Only add drawing pointer if this is the active drawing target
    if (drawingPointerPosition !== null && isDrawing && isActiveDrawingTarget) {
      allPoints.push(drawingPointerPosition);
    }
    if (allPoints.length < 2) {
      setLineSegments([]);
      setVirtualMarkerPositions([]);
      setCurveLength(0);
      // Unregister curve if we drop below 2 points
      curveRef.current = null;
      unregisterCurveRef(geometryKey);
      return;
    }
    const curve = new THREE.CatmullRomCurve3(allPoints, false, 'catmullrom', 0.5);
    curveRef.current = curve;

    // Register curve ref with store for Camera access
    registerCurveRef(geometryKey, curve);

    const curvePoints = curve.getPoints(divisions * allPoints.length);
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

  }, [markerPositions, drawingPointerPosition, isDrawing, isActiveDrawingTarget, divisions, setCurveLength, geometryKey, registerCurveRef, unregisterCurveRef]);

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
      color,
      material,
      variant: fullData.variant, // variant is not destructured, keep fullData reference
      divisions,
      thickness,
      marker,
      virtual_marker,
    };

    try {
      await createGeometry(roomId, geometryKey, "Curve", geometryData);
    } catch (error) {
      console.error("Error updating geometry:", error);
    }
  }, [roomId, markerPositions, geometryKey, color, material, variant, divisions, thickness, marker, virtual_marker, clientId]);

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

  // Early returns after all hooks (Rules of Hooks)
  if (!roomId) return null;
  if (!marker || !virtual_marker) return null;

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
          color={(isDrawing && isActiveDrawingTarget && !drawingIsValid) ? "#FF0000" : getLineColor(color)}
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
              color={(isDrawing && isActiveDrawingTarget && !drawingIsValid) ? "#FF0000" : getMarkerColor(marker.color)}
              opacity={marker.opacity}
              transparent={marker.opacity < 1}
            />
          </Dodecahedron>
        ))}

      {/* Virtual Markers */}
      {virtual_marker.enabled && virtualMarkerPositions.map((position, index) => (
          <Dodecahedron
            key={`virtual-${index}`}
            position={position}
            args={[virtual_marker.size]}
            onClick={() => handleVirtualMarkerClick(index)}
          >
            <meshBasicMaterial
              color={getMarkerColor(virtual_marker.color)}
              opacity={virtual_marker.opacity}
              transparent={virtual_marker.opacity < 1}
            />
          </Dodecahedron>
      ))}

      {/* Drawing Marker - only show for active drawing target */}
      {isDrawing && isActiveDrawingTarget && drawingPointerPosition && marker.enabled && (
        <Dodecahedron
          key={`virtual-drawing`}
          position={drawingPointerPosition}
          args={[marker.size]}
          onClick={handleDrawingMarkerClick}
        >
          <meshBasicMaterial
            color={(isDrawing && isActiveDrawingTarget && !drawingIsValid) ? "#FF0000" : getMarkerColor(marker.color)}
            opacity={marker.opacity}
            transparent={marker.opacity < 1}
          />
        </Dodecahedron>
      )}
    </group>
  );
}