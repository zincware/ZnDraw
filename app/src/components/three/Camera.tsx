import { useEffect, useMemo, useState } from "react";
import * as THREE from "three";
import { useThree } from "@react-three/fiber";
import { useAppStore } from "../../store";
import {
	isCurveAttachment,
	CurveAttachment,
	PositionType,
} from "../../utils/cameraUtils";

interface CameraData {
	// Position and target can be either direct coordinates or CurveAttachment
	position: PositionType;
	target: PositionType;

	// Camera properties
	up: [number, number, number];
	camera_type: "PerspectiveCamera" | "OrthographicCamera";
	fov: number;
	near: number;
	far: number;
	zoom: number;

	// Helper visualization
	helper_visible: boolean;
	helper_color: string;
}

/**
 * Camera geometry component for ZnDraw.
 * Renders a camera helper and handles curve attachments for position/target.
 */
export default function Camera({
	data,
	geometryKey,
}: {
	data: CameraData;
	geometryKey: string;
}) {
	const { camera: sceneCamera } = useThree();

	const { roomId, geometries, attachedCameraKey, curveRefs } = useAppStore();

	const [computedPosition, setComputedPosition] = useState<THREE.Vector3>(
		new THREE.Vector3(0, 0, 10),
	);
	const [computedTarget, setComputedTarget] = useState<THREE.Vector3>(
		new THREE.Vector3(0, 0, 0),
	);

	const isAttached = attachedCameraKey === geometryKey;

	// Extract curve info from position/target (either direct coords or CurveAttachment)
	const positionCurveKey = isCurveAttachment(data.position)
		? data.position.geometry_key
		: null;
	const positionProgress = isCurveAttachment(data.position)
		? data.position.progress
		: 0;
	const targetCurveKey = isCurveAttachment(data.target)
		? data.target.geometry_key
		: null;
	const targetProgress = isCurveAttachment(data.target)
		? data.target.progress
		: 0;

	// Get curve refs from store (shared with Curve components)
	const positionCurve = positionCurveKey
		? curveRefs[positionCurveKey]
		: undefined;
	const targetCurve = targetCurveKey ? curveRefs[targetCurveKey] : undefined;

	/**
	 * Helper: Resolve position from either direct coordinates or CurveAttachment
	 */
	const resolvePositionToVector = (
		positionData: PositionType,
		curveKey: string | null,
		curve: THREE.CatmullRomCurve3 | undefined,
		progress: number,
		fallback: THREE.Vector3,
	): THREE.Vector3 => {
		// Direct coordinates - just use them
		if (Array.isArray(positionData)) {
			return new THREE.Vector3(
				positionData[0],
				positionData[1],
				positionData[2],
			);
		}

		// CurveAttachment - resolve via curve
		if (!curveKey) {
			return fallback;
		}

		if (curve) {
			// Multi-point curve - use shared THREE.js curve object
			return curve.getPointAt(progress);
		} else {
			// Single-point curve (or not yet built) - read position directly from geometry data
			const curveGeometry = geometries[curveKey];
			if (
				curveGeometry?.type === "Curve" &&
				curveGeometry.data?.position?.[0]
			) {
				const [x, y, z] = curveGeometry.data.position[0];
				return new THREE.Vector3(x, y, z);
			} else {
				console.warn(
					`Camera ${geometryKey}: curve key '${curveKey}' not found or invalid`,
				);
				return fallback;
			}
		}
	};

	// Compute position (from direct coords or curve)
	useEffect(() => {
		const point = resolvePositionToVector(
			data.position,
			positionCurveKey,
			positionCurve,
			positionProgress,
			new THREE.Vector3(0, 0, 10),
		);
		setComputedPosition(point);
	}, [
		data.position,
		positionCurve,
		positionCurveKey,
		positionProgress,
		geometries,
		geometryKey,
	]);

	// Compute target (from direct coords or curve)
	useEffect(() => {
		const point = resolvePositionToVector(
			data.target,
			targetCurveKey,
			targetCurve,
			targetProgress,
			new THREE.Vector3(0, 0, 0),
		);
		setComputedTarget(point);
	}, [
		data.target,
		targetCurve,
		targetCurveKey,
		targetProgress,
		geometries,
		geometryKey,
	]);

	// Update scene camera if this camera is attached
	useEffect(() => {
		if (!isAttached) return;

		sceneCamera.position.copy(computedPosition);
		sceneCamera.up.set(data.up[0], data.up[1], data.up[2]);
		sceneCamera.lookAt(computedTarget);

		// Update projection properties
		let needsProjectionUpdate = false;

		if ("fov" in sceneCamera) {
			(sceneCamera as THREE.PerspectiveCamera).fov = data.fov;
			needsProjectionUpdate = true;
		}
		if ("near" in sceneCamera) {
			sceneCamera.near = data.near;
			needsProjectionUpdate = true;
		}
		if ("far" in sceneCamera) {
			sceneCamera.far = data.far;
			needsProjectionUpdate = true;
		}
		if ("zoom" in sceneCamera) {
			sceneCamera.zoom = data.zoom;
			needsProjectionUpdate = true;
		}

		// Update projection matrix if any projection property changed
		if (needsProjectionUpdate) {
			sceneCamera.updateProjectionMatrix();
		}
	}, [
		isAttached,
		computedPosition,
		computedTarget,
		data.up,
		data.fov,
		data.near,
		data.far,
		data.zoom,
		sceneCamera,
	]);

	// Create a helper camera for visualization
	const helperCamera = useMemo(() => {
		if (data.camera_type === "PerspectiveCamera") {
			return new THREE.PerspectiveCamera(data.fov, 1.0, data.near, data.far);
		} else {
			return new THREE.OrthographicCamera(-1, 1, 1, -1, data.near, data.far);
		}
	}, [data.camera_type, data.fov, data.near, data.far]);

	// Update helper camera position and target
	useEffect(() => {
		if (!helperCamera) return;

		helperCamera.position.copy(computedPosition);
		helperCamera.up.set(data.up[0], data.up[1], data.up[2]);
		helperCamera.lookAt(computedTarget);
		helperCamera.updateMatrixWorld(true);
	}, [helperCamera, computedPosition, computedTarget, data.up]);

	// Create camera helper for visualization
	const cameraHelper = useMemo(() => {
		if (helperCamera) {
			return new THREE.CameraHelper(helperCamera);
		}
		return null;
	}, [helperCamera]);

	// Update helper when camera moves
	useEffect(() => {
		if (cameraHelper) {
			cameraHelper.update();
		}
	}, [cameraHelper, computedPosition, computedTarget]);

	if (!roomId || !data.helper_visible) return null;

	return (
		<group>
			{/* Camera helper visualization (cone) */}
			{cameraHelper && <primitive object={cameraHelper} />}

			{/* Target marker (small sphere) */}
			<mesh position={computedTarget}>
				<sphereGeometry args={[0.1, 16, 16]} />
				<meshBasicMaterial
					color={data.helper_color}
					opacity={0.5}
					transparent
				/>
			</mesh>

			{/* Line from camera to target */}
			<line>
				<bufferGeometry>
					<bufferAttribute
						attach="attributes-position"
						count={2}
						itemSize={3}
						args={[
							new Float32Array([
								computedPosition.x,
								computedPosition.y,
								computedPosition.z,
								computedTarget.x,
								computedTarget.y,
								computedTarget.z,
							]),
							3,
						]}
					/>
				</bufferGeometry>
				<lineBasicMaterial
					color={data.helper_color}
					opacity={0.3}
					transparent
				/>
			</line>
		</group>
	);
}
