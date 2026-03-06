import { useColorScheme } from "@mui/material/styles";
import { ContactShadows } from "@react-three/drei";
import { useThree } from "@react-three/fiber";
import { useEffect, useRef } from "react";
import * as THREE from "three";
import { useAppStore } from "../../store";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

interface FloorData {
	active: boolean;
	position: [number, number, number];
	color: string;
	size: number;
	height: number;
	grid_spacing: number;
	grid_color: string;
	grid_opacity: number;
	show_grid: boolean;
	show_shadows: boolean;
	shadow_opacity: number;
	shadow_blur: number;
	fog_enabled: boolean;
	fog_color: string;
	fog_near: number;
	fog_far: number;
}

export const Floor = ({ data }: { data: FloorData }) => {
	const { scene } = useThree();
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);
	const gridRef = useRef<THREE.GridHelper | null>(null);
	const { mode } = useColorScheme();

	// Merge with defaults from Pydantic (single source of truth)
	const fullData = getGeometryWithDefaults<FloorData>(
		data,
		"Floor",
		geometryDefaults,
	);

	// Cap shadow scale to maintain quality regardless of floor size
	const shadowScale = Math.min(fullData.size * 0.3, 150);

	// Map "default" colors to theme colors
	const floorColor =
		fullData.color === "default"
			? mode === "light"
				? "#e0e0e0"
				: "#303030"
			: fullData.color;

	const gridColor =
		fullData.grid_color === "default"
			? mode === "light"
				? "#616161"
				: "#9e9e9e"
			: fullData.grid_color;

	const fogColor =
		fullData.fog_color === "default"
			? mode === "light"
				? "#FFFFFF"
				: "#212121"
			: fullData.fog_color;

	// Sync fog with floor settings
	useEffect(() => {
		if (!fullData.fog_enabled) {
			scene.fog = null;
			return;
		}

		scene.fog = new THREE.Fog(fogColor, fullData.fog_near, fullData.fog_far);

		return () => {
			scene.fog = null;
		};
	}, [
		fullData.fog_enabled,
		fullData.fog_near,
		fullData.fog_far,
		fogColor,
		scene,
	]);

	// Update grid when settings change
	useEffect(() => {
		if (!gridRef.current || !fullData.show_grid) return;

		const grid = gridRef.current;
		grid.position.y = fullData.height;

		// Update grid material color/opacity
		if (grid.material instanceof THREE.Material) {
			grid.material.opacity = fullData.grid_opacity;
			grid.material.transparent = true;
			grid.material.color.set(gridColor);
		}
	}, [fullData.height, gridColor, fullData.grid_opacity, fullData.show_grid]);

	const divisions = Math.floor(fullData.size / fullData.grid_spacing);

	return (
		<group>
			{/* Floor Plane (receives shadows) */}
			<mesh
				rotation={[-Math.PI / 2, 0, 0]}
				position={[0, fullData.height - 0.01, 0]}
				receiveShadow
			>
				<planeGeometry args={[fullData.size, fullData.size]} />
				<meshStandardMaterial
					color={floorColor}
					roughness={0.8}
					metalness={0.2}
				/>
			</mesh>

			{/* Grid Helper */}
			{fullData.show_grid && (
				<gridHelper
					ref={gridRef}
					args={[fullData.size, divisions, gridColor, gridColor]}
					position={[0, fullData.height, 0]}
				/>
			)}

			{/* Contact Shadows (drei component) */}
			{fullData.show_shadows && (
				<ContactShadows
					position={[0, fullData.height + 0.01, 0]}
					opacity={fullData.shadow_opacity}
					scale={shadowScale}
					blur={fullData.shadow_blur}
					far={20}
					resolution={512}
				/>
			)}
		</group>
	);
};
