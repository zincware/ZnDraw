import { useTheme } from "@mui/material/styles";
import { Line } from "@react-three/drei";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import { useGeometryEditing } from "../../hooks/useGeometryEditing";
import { useGeometryPersistence } from "../../hooks/useGeometryPersistence";
import { useAppStore } from "../../store";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

interface InteractionSettings {
	enabled: boolean;
	color: string;
	opacity: number;
}

interface CircleCurveData {
	position: [number, number, number][];
	radius: number;
	start_angle: number;
	end_angle: number;
	rotation: [number, number, number][];
	scale: [number, number, number][];
	color: string;
	material: "LineBasicMaterial" | "LineDashedMaterial";
	divisions: number;
	thickness: number;
	selecting?: InteractionSettings;
	hovering?: InteractionSettings;
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
	const updateSelections = useAppStore((state) => state.updateSelections);
	const selections = useAppStore((state) => state.selections);
	const hoveredGeometryInstance = useAppStore(
		(state) => state.hoveredGeometryInstance,
	);
	const setHoveredGeometryInstance = useAppStore(
		(state) => state.setHoveredGeometryInstance,
	);

	const theme = useTheme();

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
		selecting,
		hovering,
	} = fullData;

	// Extract single values from length-1 arrays
	const pos = position[0];
	const rot = rotation[0];
	const scl = scale[0];

	const groupRef = useRef<THREE.Group>(null);
	const [isHovered, setIsHovered] = useState(false);

	// Resolve "default" color via theme
	const resolvedColor =
		color === "default" ? theme.palette.text.primary : color;

	// Selection state for this geometry (single instance, index 0)
	const selectedIndices = useMemo(
		() => (selections[geometryKey] || []).filter((i: number) => i === 0),
		[selections, geometryKey],
	);
	const isSelected = selectedIndices.length > 0;

	// Hover state
	const isHoveredInstance =
		hoveredGeometryInstance?.geometryKey === geometryKey &&
		hoveredGeometryInstance?.instanceId === 0;

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
		pos,
		rot,
		scl,
		geometryKey,
		registerCurveRef,
		unregisterCurveRef,
	]);

	// Persistence: save geometry changes back to server
	useGeometryPersistence(geometryKey, "CircleCurve");

	// Editing: plug into transform controls system
	useGeometryEditing(
		geometryKey,
		position,
		rotation,
		scale,
		selectedIndices,
		"CircleCurve",
		fullData,
		1,
	);

	// Marker event handlers
	const handleClick = useCallback(
		(event: any) => {
			if (event.detail !== 1) return;
			event.stopPropagation();
			updateSelections(geometryKey, 0, event.shiftKey);
		},
		[geometryKey, updateSelections],
	);

	const handlePointerOver = useCallback(
		(event: any) => {
			event.stopPropagation();
			setIsHovered(true);
			setHoveredGeometryInstance(geometryKey, 0);
		},
		[geometryKey, setHoveredGeometryInstance],
	);

	const handlePointerOut = useCallback(() => {
		setIsHovered(false);
		setHoveredGeometryInstance(null, null);
	}, [setHoveredGeometryInstance]);

	// Determine marker color based on state
	const markerColor = useMemo(() => {
		if (isSelected && selecting?.enabled) return selecting.color;
		if ((isHovered || isHoveredInstance) && hovering?.enabled)
			return hovering.color;
		return theme.palette.secondary.main;
	}, [
		isSelected,
		isHovered,
		isHoveredInstance,
		selecting,
		hovering,
		theme.palette.secondary.main,
	]);

	const markerOpacity = useMemo(() => {
		if (isSelected && selecting?.enabled) return selecting.opacity;
		if ((isHovered || isHoveredInstance) && hovering?.enabled)
			return hovering.opacity;
		return 0.5;
	}, [isSelected, isHovered, isHoveredInstance, selecting, hovering]);

	// Marker geometry (shared)
	const markerGeometry = useMemo(() => new THREE.DodecahedronGeometry(1), []);

	// Hide line visuals during pathtracing or when inactive
	// But keep mounted so curve ref stays registered for CurveAttachment
	if (pathtracingEnabled || fullData.active === false) {
		return (
			<group
				ref={groupRef}
				position={pos}
				rotation={new THREE.Euler(...rot)}
				scale={scl}
			/>
		);
	}

	if (points.length < 2) return null;

	return (
		<group
			ref={groupRef}
			position={pos}
			rotation={new THREE.Euler(...rot)}
			scale={scl}
		>
			<Line
				points={points}
				color={resolvedColor}
				lineWidth={thickness}
				dashed={material === "LineDashedMaterial"}
			/>

			{/* Center marker for click/hover/selection */}
			<mesh
				geometry={markerGeometry}
				scale={0.5}
				onClick={handleClick}
				onPointerOver={handlePointerOver}
				onPointerOut={handlePointerOut}
			>
				<meshBasicMaterial
					color={markerColor}
					opacity={markerOpacity}
					transparent={markerOpacity < 1}
				/>
			</mesh>
		</group>
	);
}
