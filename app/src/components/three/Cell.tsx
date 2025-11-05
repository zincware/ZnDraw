import { useColorScheme, useTheme } from "@mui/material/styles";
import { Line } from "@react-three/drei";
import { useAppStore } from "../../store";
import { useState, useEffect } from "react";
import * as THREE from "three";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { getFrames } from "../../myapi/client";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

interface CellData {
  position: string;
  color: string;
  material: "LineBasicMaterial" | "LineDashedMaterial";
  thickness: number;
}

export const Cell = ({ data }: {data: CellData}) => {
  const { mode } = useColorScheme();
  const theme = useTheme();
  // Use individual selectors to prevent unnecessary re-renders
  const currentFrame = useAppStore((state) => state.currentFrame);
  const roomId = useAppStore((state) => state.roomId);
  const userName = useAppStore((state) => state.userName);
  const geometryDefaults = useAppStore((state) => state.geometryDefaults);

  // Merge with defaults from Pydantic (single source of truth)
  const fullData = getGeometryWithDefaults<CellData>(data, "Cell", geometryDefaults);

  const [displayedVertices, setDisplayedVertices] = useState<
    THREE.Vector3[][] | null
  >(null);

  // Simple single query for cell data
  const { data: cellData, isFetching } = useQuery({
    queryKey: ["frame", roomId, currentFrame, fullData.position],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [fullData.position], signal),
    enabled: !!roomId && !!userName,
    placeholderData: keepPreviousData,
  });

  useEffect(() => {
    if (isFetching) {
      return; // Wait for data, keepPreviousData ensures old view remains
    }

    const cell = cellData?.[fullData.position];
    if (!cell) {
      setDisplayedVertices(null);
      return;
    }

    // cell is directly the TypedArray (Float32Array or Float64Array)
    // Expected: 9 values representing a 3x3 matrix
    if (cell.length !== 9) {
      console.error("Invalid cell data length:", cell.length, "expected 9");
      setDisplayedVertices(null);
      return;
    }

    const cellVectors = [
      new THREE.Vector3(cell[0], cell[1], cell[2]),
      new THREE.Vector3(cell[3], cell[4], cell[5]),
      new THREE.Vector3(cell[6], cell[7], cell[8]),
    ];

    const origin = new THREE.Vector3(0, 0, 0);
    const v = [
      origin,
      cellVectors[0],
      cellVectors[1],
      cellVectors[0].clone().add(cellVectors[1]),
      cellVectors[2],
      cellVectors[0].clone().add(cellVectors[2]),
      cellVectors[1].clone().add(cellVectors[2]),
      cellVectors[0].clone().add(cellVectors[1]).add(cellVectors[2]),
    ];

    const calculatedVertices = [
      [v[0], v[1]],
      [v[1], v[3]],
      [v[3], v[2]],
      [v[2], v[0]],
      [v[4], v[5]],
      [v[5], v[7]],
      [v[7], v[6]],
      [v[6], v[4]],
      [v[0], v[4]],
      [v[1], v[5]],
      [v[2], v[6]],
      [v[3], v[7]],
    ];
    setDisplayedVertices(calculatedVertices);
  }, [cellData, isFetching, fullData.position]);

  const lineColor = fullData.color === "default" ? (mode === "light" ? theme.palette.primary.dark : theme.palette.primary.light) : fullData.color;

  return (
    <group>
      {displayedVertices &&
        displayedVertices.map((points, index) => (
          <Line key={index} points={points} color={lineColor} lineWidth={fullData.thickness} dashed={fullData.material === "LineDashedMaterial"}/>
        ))}
    </group>
  );
};
