import { useColorScheme, useTheme } from "@mui/material/styles";
import { Line } from "@react-three/drei";
import {
	keepPreviousData,
	useQuery,
	useQueryClient,
} from "@tanstack/react-query";
import { useEffect, useMemo, useState } from "react";
import * as THREE from "three";
import { getFrameBatched } from "../../hooks/useFrameBatch";
import { useAppStore } from "../../store";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

interface CellData {
	position: string;
	color: string;
	material: "LineBasicMaterial" | "LineDashedMaterial";
	thickness: number;
}

/**
 * Renders the periodic cell as line segments.
 *
 * Note: Cell uses Line components which are not supported by the GPU pathtracer.
 * When pathtracingEnabled is true, the cell is hidden.
 */
export const Cell = ({
	data,
	pathtracingEnabled = false,
}: { data: CellData; pathtracingEnabled?: boolean }) => {
	const queryClient = useQueryClient();
	const { mode } = useColorScheme();
	const theme = useTheme();
	// Use individual selectors to prevent unnecessary re-renders
	const currentFrame = useAppStore((state) => state.currentFrame);
	const frameCount = useAppStore((state) => state.frameCount);
	const roomId = useAppStore((state) => state.roomId);
	const userName = useAppStore((state) => state.user?.email ?? null);
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);

	// Merge with defaults from Pydantic (single source of truth)
	const fullData = useMemo(
		() => getGeometryWithDefaults<CellData>(data, "Cell", geometryDefaults),
		[data, geometryDefaults],
	);

	const [displayedVertices, setDisplayedVertices] = useState<
		THREE.Vector3[][] | null
	>(null);

	// Simple single query for cell data
	const { data: cellData, isFetching } = useQuery({
		queryKey: ["frame", roomId, currentFrame, fullData.position],
		queryFn: () =>
			getFrameBatched(roomId!, currentFrame, fullData.position, queryClient),
		enabled: !!roomId && !!userName,
		placeholderData: keepPreviousData,
	});

	useEffect(() => {
		// When frameCount is 0, explicitly clear cell (e.g., after del vis[:])
		if (frameCount === 0) {
			if (displayedVertices !== null) setDisplayedVertices(null);
			return;
		}

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
			new THREE.Vector3(Number(cell[0]), Number(cell[1]), Number(cell[2])),
			new THREE.Vector3(Number(cell[3]), Number(cell[4]), Number(cell[5])),
			new THREE.Vector3(Number(cell[6]), Number(cell[7]), Number(cell[8])),
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
	}, [cellData, frameCount, isFetching, fullData.position]);

	const lineColor =
		fullData.color === "default"
			? mode === "light"
				? theme.palette.primary.dark
				: theme.palette.primary.light
			: fullData.color;

	// Hide cell when pathtracing (Line components not supported by GPU pathtracer)
	if (pathtracingEnabled) return null;

	return (
		<group>
			{displayedVertices &&
				displayedVertices.map((points, index) => (
					<Line
						key={index}
						points={points}
						color={lineColor}
						lineWidth={fullData.thickness}
						dashed={fullData.material === "LineDashedMaterial"}
					/>
				))}
		</group>
	);
};
