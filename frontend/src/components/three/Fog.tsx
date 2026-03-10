/**
 * Fog component for distance-based atmospheric fog effect.
 *
 * Fog is a scene object (geometry) that fades objects to a color based on distance.
 * Use 'default' color to automatically match the theme background.
 */

import { useColorScheme } from "@mui/material/styles";
import { useThree } from "@react-three/fiber";
import { useEffect } from "react";
import * as THREE from "three";
import { useAppStore } from "../../store";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

interface FogData {
	active: boolean;
	color: string;
	near: number;
	far: number;
}

interface FogProps {
	geometryKey: string;
	data: FogData;
}

export function Fog({ data }: FogProps) {
	const { scene } = useThree();
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);
	const { mode } = useColorScheme();

	const fullData = getGeometryWithDefaults<FogData>(
		data,
		"Fog",
		geometryDefaults,
	);

	// Map "default" color to theme background
	const fogColor =
		fullData.color === "default"
			? mode === "light"
				? "#FFFFFF"
				: "#212121"
			: fullData.color;

	// Sync fog with scene
	useEffect(() => {
		if (!fullData.active) {
			scene.fog = null;
			return;
		}

		scene.fog = new THREE.Fog(fogColor, fullData.near, fullData.far);

		return () => {
			scene.fog = null;
		};
	}, [fullData.active, fullData.near, fullData.far, fogColor, scene]);

	// This component doesn't render anything visible - it just manages scene.fog
	return null;
}
