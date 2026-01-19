import { useEffect, useMemo, useState } from "react";
import * as THREE from "three";
import { useThree } from "@react-three/fiber";
import { useAppStore } from "../../store";
import {
	isCurveAttachment,
	resolvePosition,
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

	const roomId = useAppStore((state) => state.roomId);
	const geometries = useAppStore((state) => state.geometries);
	const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);
	const curveRefs = useAppStore((state) => state.curveRefs);

	// Initial position/target from shared utility (linear interpolation)
	const [computedPosition, setComputedPosition] = useState<THREE.Vector3>(
		() => {
			const [x, y, z] = resolvePosition(data.position, geometries);
			return new THREE.Vector3(x, y, z);
		},
	);
	const [computedTarget, setComputedTarget] = useState<THREE.Vector3>(() => {
		const [x, y, z] = resolvePosition(data.target, geometries);
		return new THREE.Vector3(x, y, z);
	});

	const isAttached = attachedCameraKey === geometryKey;

	// Extract curve info from position/target
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
	 * Resolve position to Vector3. Uses curveRef for smooth spline interpolation
	 * when available, falls back to resolvePosition (linear) otherwise.
	 */
	const resolveToVector3 = (
		positionData: PositionType,
		curve: THREE.CatmullRomCurve3 | undefined,
		progress: number,
	): THREE.Vector3 => {
		// Use curveRef for smooth spline interpolation when available
		if (curve && isCurveAttachment(positionData)) {
			return curve.getPointAt(progress);
		}
		// Fallback to linear interpolation via shared utility
		const [x, y, z] = resolvePosition(positionData, geometries);
		return new THREE.Vector3(x, y, z);
	};

	// Compute position (from direct coords or curve)
	useEffect(() => {
		const point = resolveToVector3(
			data.position,
			positionCurve,
			positionProgress,
		);
		setComputedPosition(point);
	}, [data.position, positionCurve, positionProgress, geometries]);

	// Compute target (from direct coords or curve)
	useEffect(() => {
		const point = resolveToVector3(data.target, targetCurve, targetProgress);
		setComputedTarget(point);
	}, [data.target, targetCurve, targetProgress, geometries]);

	// Update scene camera if this camera is attached
	useEffect(() => {
		if (!isAttached) return;

		sceneCamera.position.copy(computedPosition);
		sceneCamera.up.set(data.up[0], data.up[1], data.up[2]);
		sceneCamera.lookAt(computedTarget);

		// Update projection properties
		if ("fov" in sceneCamera) {
			(sceneCamera as THREE.PerspectiveCamera).fov = data.fov;
		}
		sceneCamera.near = data.near;
		sceneCamera.far = data.far;
		sceneCamera.zoom = data.zoom;
		sceneCamera.updateProjectionMatrix();
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
		}
		return new THREE.OrthographicCamera(-1, 1, 1, -1, data.near, data.far);
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
