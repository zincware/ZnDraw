import { Line } from "@react-three/drei";
import { useEffect, useMemo, useRef } from "react";
import * as THREE from "three";
import { useGeometryPersistence } from "../../hooks/useGeometryPersistence";
import { useAppStore } from "../../store";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

interface CircleCurveData {
	position: [number, number, number];
	radius: number;
	start_angle: number;
	end_angle: number;
	rotation: [number, number, number];
	scale: [number, number, number];
	color: string;
	material: "LineBasicMaterial" | "LineDashedMaterial";
	divisions: number;
	thickness: number;
	active?: boolean;
}

/**
 * A custom THREE.Curve that produces 3D points on an ellipse in the XY plane.
 * Used so CurveAttachment can call getPointAt(t) to get world-space positions.
 */
class EllipseCurve3D extends THREE.Curve<THREE.Vector3> {
	constructor(
		private radius: number,
		private startRad: number,
		private endRad: number,
		private groupRef: React.RefObject<THREE.Group | null>,
	) {
		super();
	}

	getPoint(t: number, optionalTarget = new THREE.Vector3()): THREE.Vector3 {
		const angle = this.startRad + t * (this.endRad - this.startRad);
		// Local point in XY plane
		const localPoint = optionalTarget.set(
			this.radius * Math.cos(angle),
			this.radius * Math.sin(angle),
			0,
		);
		// Transform to world space via the group's world matrix
		if (this.groupRef.current) {
			this.groupRef.current.updateWorldMatrix(true, false);
			localPoint.applyMatrix4(this.groupRef.current.matrixWorld);
		}
		return localPoint;
	}
}

export default function CircleCurve({
	data,
	geometryKey,
	pathtracingEnabled = false,
}: {
	data: CircleCurveData;
	geometryKey: string;
	pathtracingEnabled?: boolean;
}) {
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);
	const registerCurveRef = useAppStore((state) => state.registerCurveRef);
	const unregisterCurveRef = useAppStore((state) => state.unregisterCurveRef);

	const fullData = useMemo(
		() =>
			getGeometryWithDefaults<CircleCurveData>(
				data,
				"CircleCurve",
				geometryDefaults,
			),
		[data, geometryDefaults],
	);

	const {
		position,
		radius,
		start_angle,
		end_angle,
		rotation,
		scale,
		color,
		material,
		divisions,
		thickness,
	} = fullData;

	const groupRef = useRef<THREE.Group>(null);

	// Convert percentage angles to radians
	const startRad = (start_angle / 100) * Math.PI * 2;
	const endRad = (end_angle / 100) * Math.PI * 2;

	// Generate 2D ellipse points and convert to 3D (local space)
	const points = useMemo(() => {
		const curve = new THREE.EllipseCurve(
			0,
			0, // center (local)
			radius,
			radius, // radii (scale handles ellipse)
			startRad,
			endRad,
			false, // not clockwise
			0, // no rotation (handled by group)
		);
		const pts2D = curve.getPoints(divisions);
		return pts2D.map((p) => new THREE.Vector3(p.x, p.y, 0));
	}, [radius, startRad, endRad, divisions]);

	// Create and register 3D curve for CurveAttachment
	useEffect(() => {
		const curve3D = new EllipseCurve3D(radius, startRad, endRad, groupRef);
		registerCurveRef(geometryKey, curve3D);
		return () => unregisterCurveRef(geometryKey);
	}, [
		radius,
		startRad,
		endRad,
		position,
		rotation,
		scale,
		geometryKey,
		registerCurveRef,
		unregisterCurveRef,
	]);

	// Persistence: save geometry changes back to server
	useGeometryPersistence(geometryKey, "CircleCurve");

	// Hide line visuals during pathtracing (Line not supported by GPU pathtracer)
	// But keep mounted so curve ref stays registered for CurveAttachment
	if (pathtracingEnabled) return <group />;

	if (points.length < 2) return null;

	return (
		<group
			ref={groupRef}
			position={position}
			rotation={new THREE.Euler(...rotation)}
			scale={scale}
		>
			<Line
				points={points}
				color={color}
				lineWidth={thickness}
				dashed={material === "LineDashedMaterial"}
			/>
		</group>
	);
}
