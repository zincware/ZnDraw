/**
 * Isosurface component for volumetric data visualization.
 *
 * Fetches mesh data from a dedicated server endpoint that runs marching cubes.
 * Does NOT use useRegisterFrameKeys or getFrameBatched — volumetric data stays
 * server-side. Only the extracted mesh (kilobytes) crosses the wire.
 */

import { useQuery } from "@tanstack/react-query";
import { useEffect, useMemo, useRef } from "react";
import * as THREE from "three";
import { useAppStore } from "../../store";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";
import { fetchIsosurface } from "../../myapi/client";

interface IsosurfaceData {
	active: boolean;
	cube_key: string;
	isovalue: number;
	color: string;
	resolution: number; // 0..1 float, converted to step_size server-side
	opacity: number;
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

	const geometryRef = useRef<THREE.BufferGeometry | null>(null);

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

	useEffect(() => {
		const prev = geometryRef.current;
		geometryRef.current = geometry;
		if (prev && prev !== geometry) {
			prev.dispose();
		}
	}, [geometry]);

	useEffect(() => {
		return () => {
			geometryRef.current?.dispose();
		};
	}, []);

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
