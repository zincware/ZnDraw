/**
 * Light components for scene illumination.
 *
 * Lights are scene objects (geometries) that can be attached to the camera
 * to follow view direction, providing consistent lighting regardless of view angle.
 */

import { useFrame, useThree } from "@react-three/fiber";
import { useRef } from "react";
import type * as THREE from "three";
import { useAppStore } from "../../store";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

// ==================== Type Definitions ====================

interface LightPosition {
	camera_attached: boolean;
	x: number;
	y: number;
	z: number;
}

interface DirectionalLightData {
	active: boolean;
	intensity: number;
	position: LightPosition;
	color: string;
}

interface AmbientLightData {
	active: boolean;
	intensity: number;
	color: string;
}

interface HemisphereLightData {
	active: boolean;
	intensity: number;
	sky_color: string;
	ground_color: string;
}

// ==================== DirectionalLight ====================

interface DirectionalLightProps {
	geometryKey: string;
	data: DirectionalLightData;
}

export function DirectionalLight({ data }: DirectionalLightProps) {
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);
	const rigRef = useRef<THREE.Group>(null!);
	const { camera } = useThree();

	const fullData = getGeometryWithDefaults<DirectionalLightData>(
		data,
		"DirectionalLight",
		geometryDefaults,
	);

	const { position } = fullData;
	const coords: [number, number, number] = [position.x, position.y, position.z];

	// If attached to camera, follow camera position/rotation each frame
	useFrame(() => {
		if (rigRef.current && position.camera_attached) {
			rigRef.current.position.copy(camera.position);
			rigRef.current.quaternion.copy(camera.quaternion);
		}
	});

	if (!fullData.active) return null;

	// Camera-attached light: wrap in group that follows camera
	if (position.camera_attached) {
		return (
			<group ref={rigRef}>
				<directionalLight
					position={coords}
					intensity={fullData.intensity}
					color={fullData.color}
				/>
			</group>
		);
	}

	// World-space light: render at fixed position
	return (
		<directionalLight
			position={coords}
			intensity={fullData.intensity}
			color={fullData.color}
		/>
	);
}

// ==================== AmbientLight ====================

interface AmbientLightProps {
	geometryKey: string;
	data: AmbientLightData;
}

export function AmbientLight({ data }: AmbientLightProps) {
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);

	const fullData = getGeometryWithDefaults<AmbientLightData>(
		data,
		"AmbientLight",
		geometryDefaults,
	);

	if (!fullData.active) return null;

	return <ambientLight intensity={fullData.intensity} color={fullData.color} />;
}

// ==================== HemisphereLight ====================

interface HemisphereLightProps {
	geometryKey: string;
	data: HemisphereLightData;
}

export function HemisphereLight({ data }: HemisphereLightProps) {
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);

	const fullData = getGeometryWithDefaults<HemisphereLightData>(
		data,
		"HemisphereLight",
		geometryDefaults,
	);

	if (!fullData.active) return null;

	return (
		<hemisphereLight
			color={fullData.sky_color}
			groundColor={fullData.ground_color}
			intensity={fullData.intensity}
		/>
	);
}
