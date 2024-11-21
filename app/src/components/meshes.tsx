import type React from "react";
import { useEffect, useMemo, useRef } from "react";
import * as THREE from "three";
import { BufferGeometryUtils } from "three/examples/jsm/Addons.js";
import { type ColorRange, type HSLColor, interpolateColor } from "./utils";
import { useMergedMesh } from "./utils/mergeInstancedMesh";

function createArrowMesh() {
	const cylinderRadius = 0.04;
	const cylinderHeight = 0.6;
	const coneRadius = 0.1;
	const coneHeight = 0.4;

	const cylinderGeometry = new THREE.CylinderGeometry(
		cylinderRadius,
		cylinderRadius,
		cylinderHeight,
		32,
	);
	const coneGeometry = new THREE.ConeGeometry(coneRadius, coneHeight, 32);

	cylinderGeometry.translate(0, cylinderHeight / 2, 0);
	coneGeometry.translate(0, cylinderHeight + coneHeight / 2, 0);

	const arrowGeometry = BufferGeometryUtils.mergeGeometries([
		cylinderGeometry,
		coneGeometry,
	]);

	return arrowGeometry;
}

interface ArrowsProps {
	start: number[][];
	end: number[][];
	scale_vector_thickness?: boolean;
	colormap: HSLColor[];
	colorrange: ColorRange;
	opacity?: number;
	rescale?: number;
	pathTracingSettings: any | undefined;
}

const Arrows: React.FC<ArrowsProps> = ({
	start,
	end,
	scale_vector_thickness,
	colormap,
	colorrange,
	opacity = 1.0,
	rescale = 1.0,
	pathTracingSettings = undefined,
}) => {
	const meshRef = useRef<THREE.InstancedMesh>(null);
	const materialRef = useRef<THREE.MeshStandardMaterial>(null);

	const geometry = useMemo(() => {
		const _geom = createArrowMesh();
		if (pathTracingSettings?.enabled) {
			// make invisible when path tracing is enabled
			_geom.scale(0, 0, 0);
		}
		return _geom;
	}, [pathTracingSettings]);

	const instancedGeometry = useMemo(() => {
		return createArrowMesh();
	}, []);

	const mergedMesh = useMergedMesh(
		meshRef,
		instancedGeometry,
		pathTracingSettings,
		[
			start,
			end,
			scale_vector_thickness,
			colormap,
			colorrange,
			opacity,
			rescale,
		],
	);

	useEffect(() => {
		if (!meshRef.current) return;
		const matrix = new THREE.Matrix4();
		const up = new THREE.Vector3(0, 1, 0);
		const startVector = new THREE.Vector3();
		const endVector = new THREE.Vector3();
		const direction = new THREE.Vector3();
		const quaternion = new THREE.Quaternion();

		for (let i = 0; i < start.length; i++) {
			startVector.fromArray(start[i]);
			endVector.fromArray(end[i]);
			direction.subVectors(endVector, startVector);
			let length = direction.length();
			const color = interpolateColor(colormap, colorrange, length);
			// rescale after the color interpolation
			length *= rescale;

			const scale = scale_vector_thickness
				? new THREE.Vector3(length, length, length)
				: new THREE.Vector3(1, length, 1);

			quaternion.setFromUnitVectors(up, direction.clone().normalize());
			matrix.makeRotationFromQuaternion(quaternion);
			matrix.setPosition(startVector);
			matrix.scale(scale);

			meshRef.current.setColorAt(i, color);
			meshRef.current.setMatrixAt(i, matrix);
		}
		meshRef.current.instanceMatrix.needsUpdate = true;
	}, [start, end, scale_vector_thickness, colormap, colorrange]);

	useEffect(() => {
		if (!materialRef.current) return;
		materialRef.current.needsUpdate = true; // TODO: check for particles as well
		if (!meshRef.current) return;
		if (!meshRef.current.instanceColor) return;
		meshRef.current.instanceColor.needsUpdate = true;
	}, [start, end, scale_vector_thickness, colormap, colorrange]);

	return (
		<instancedMesh ref={meshRef} args={[geometry, undefined, start.length]}>
			<meshStandardMaterial
				ref={materialRef}
				attach="material"
				transparent
				opacity={opacity}
			/>
		</instancedMesh>
	);
};

export default Arrows;
