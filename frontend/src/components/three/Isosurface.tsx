/**
 * Isosurface component for volumetric data visualization.
 *
 * Fetches mesh data from a dedicated server endpoint that runs marching cubes.
 * Does NOT use useRegisterFrameKeys or getFrameBatched — volumetric data stays
 * server-side. Only the extracted mesh (kilobytes) crosses the wire.
 */

import { useQuery } from "@tanstack/react-query";
import { useMemo } from "react";
import * as THREE from "three";
import { useAppStore } from "../../store";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";
import { unpackBinary } from "../../utils/msgpack-numpy";
import axios from "axios";
import { getToken } from "../../utils/auth";

interface IsosurfaceData {
	active: boolean;
	cube_key: string;
	isovalue: number;
	color: string;
	resolution: number; // 0..1 float, converted to step_size server-side
	opacity: number;
}

async function fetchIsosurface(
	roomId: string,
	frame: number,
	cubeKey: string,
	isovalue: number,
	resolution: number,
): Promise<{
	vertices: number[][] | Float32Array;
	faces: number[][] | Uint32Array;
}> {
	const params = new URLSearchParams({
		cube_key: cubeKey,
		isovalue: isovalue.toString(),
		resolution: resolution.toString(),
	});
	const token = getToken();
	const response = await axios.get(
		`/v1/rooms/${roomId}/frames/${frame}/isosurface?${params.toString()}`,
		{
			responseType: "arraybuffer",
			headers: token ? { Authorization: `Bearer ${token}` } : {},
		},
	);
	return unpackBinary(new Uint8Array(response.data)) as {
		vertices: number[][] | Float32Array;
		faces: number[][] | Uint32Array;
	};
}

export default function Isosurface({
	data,
}: {
	data: IsosurfaceData;
	geometryKey: string;
}) {
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);
	const roomId = useAppStore((state) => state.roomId);
	const currentFrame = useAppStore((state) => state.currentFrame);

	const fullData = getGeometryWithDefaults<IsosurfaceData>(
		data,
		"Isosurface",
		geometryDefaults,
	);

	const { data: meshData } = useQuery({
		queryKey: [
			"isosurface",
			roomId,
			currentFrame,
			fullData.cube_key,
			fullData.isovalue,
			fullData.resolution,
		],
		queryFn: () =>
			fetchIsosurface(
				roomId!,
				currentFrame,
				fullData.cube_key,
				fullData.isovalue,
				fullData.resolution,
			),
		enabled: !!roomId && !!fullData.cube_key && fullData.active,
		staleTime: 30000,
	});

	const geometry = useMemo(() => {
		if (!meshData || meshData.vertices.length === 0) return null;

		const geo = new THREE.BufferGeometry();

		// msgpack-numpy may return Float32Array or nested number[][]
		const vertices =
			meshData.vertices instanceof Float32Array
				? meshData.vertices
				: new Float32Array(
						Array.isArray(meshData.vertices[0])
							? (meshData.vertices as number[][]).flat()
							: (meshData.vertices as unknown as number[]),
					);
		geo.setAttribute("position", new THREE.BufferAttribute(vertices, 3));

		const indices =
			meshData.faces instanceof Uint32Array
				? meshData.faces
				: new Uint32Array(
						Array.isArray(meshData.faces[0])
							? (meshData.faces as number[][]).flat()
							: (meshData.faces as unknown as number[]),
					);
		geo.setIndex(new THREE.BufferAttribute(indices, 1));

		geo.computeVertexNormals();
		return geo;
	}, [meshData]);

	if (!fullData.active || !geometry) return null;

	const isTransparent = fullData.opacity < 1.0;

	return (
		<mesh geometry={geometry}>
			<meshPhysicalMaterial
				color={fullData.color}
				side={THREE.DoubleSide}
				transparent={isTransparent}
				opacity={fullData.opacity}
				depthWrite={!isTransparent}
				roughness={0.4}
			/>
		</mesh>
	);
}
