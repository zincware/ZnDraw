import { useColorScheme, useTheme } from "@mui/material/styles";
import { Line } from "@react-three/drei";
import { useAppStore } from "../../store";
import { useState, useEffect, useMemo } from "react";
import * as THREE from "three";
import { useQueries } from "@tanstack/react-query";
import { getFrameDataOptions } from "../../hooks/useTrajectoryData";

export const SimulationCell = () => {
  const { mode } = useColorScheme();
  const theme = useTheme();
  const { currentFrame, roomId } = useAppStore();

  // It's slightly better practice to initialize state with null instead of undefined
  const [displayedVertices, setDisplayedVertices] = useState<
    THREE.Vector3[][] | null
  >(null);

  const queries = useMemo(() => {
    if (!roomId) return [];
    return [getFrameDataOptions(roomId, currentFrame, "cell")];
  }, [roomId, currentFrame]);

  const [cellQueryResult] = useQueries({ queries });

  // 1. Destructure the values you need from the query result.
  // This makes the dependency array cleaner.
  const { data, isSuccess, isFetching } = cellQueryResult || {};

  // 2. Create a stable key from the data that only changes when the content changes.
  // This is the key to preventing unnecessary effect runs.
  const stableDataKey = useMemo(() => JSON.stringify(data), [data]);

  useEffect(() => {
    // We can still log here, but it will no longer cause a loop.
    console.log("Running effect because data changed. New key:", stableDataKey);

    if (isSuccess && data) {
      const { data: flatData, shape } = data;
      if (!shape || shape.length !== 2 || shape[0] !== 3 || shape[1] !== 3) {
        console.error("Invalid cell dimensions received from query:", shape);
        return;
      }

      const cellVectors = [
        new THREE.Vector3(flatData[0], flatData[1], flatData[2]),
        new THREE.Vector3(flatData[3], flatData[4], flatData[5]),
        new THREE.Vector3(flatData[6], flatData[7], flatData[8]),
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
    }
  }, [isSuccess, stableDataKey]); // NO LONGER DEPENDS ON `cellQueryResult`

  const lineColor =
    mode === "light" ? theme.palette.primary.dark : theme.palette.primary.light;

  return (
    <group>
      {displayedVertices &&
        displayedVertices.map((points, index) => (
          <Line key={index} points={points} color={lineColor} lineWidth={2} />
        ))}
    </group>
  );
};
