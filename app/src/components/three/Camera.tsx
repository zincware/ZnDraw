import { useEffect, useMemo, useState } from "react";
import * as THREE from "three";
import { useThree } from "@react-three/fiber";
import { useAppStore } from "../../store";

interface CameraData {
  // Curve references (always used)
  position_curve_key: string;
  position_progress: number;
  target_curve_key: string;
  target_progress: number;

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

  const {
    roomId,
    geometries,
    attachedCameraKey,
    curveRefs,
  } = useAppStore();

  const [computedPosition, setComputedPosition] = useState<THREE.Vector3>(
    new THREE.Vector3(0, 0, 10)
  );
  const [computedTarget, setComputedTarget] = useState<THREE.Vector3>(
    new THREE.Vector3(0, 0, 0)
  );

  const isAttached = attachedCameraKey === geometryKey;

  // Get curve refs from store (shared with Curve components)
  const positionCurveKey = data.position_curve_key || null;
  const targetCurveKey = data.target_curve_key || null;
  const positionCurve = positionCurveKey ? curveRefs[positionCurveKey] : undefined;
  const targetCurve = targetCurveKey ? curveRefs[targetCurveKey] : undefined;

  /**
   * Helper: Compute 3D point from curve key + progress
   * Handles both multi-point curves (via THREE.js curve object) and single-point curves (direct read)
   */
  const computePointFromCurve = (
    curveKey: string | null,
    curve: THREE.CatmullRomCurve3 | undefined,
    progress: number,
    fallback: THREE.Vector3
  ): THREE.Vector3 => {
    if (!curveKey) {
      return fallback;
    }

    if (curve) {
      // Multi-point curve - use shared THREE.js curve object
      return curve.getPointAt(progress);
    } else {
      // Single-point curve (or not yet built) - read position directly from geometry data
      const curveGeometry = geometries[curveKey];
      if (curveGeometry?.type === "Curve" && curveGeometry.data?.position?.[0]) {
        const [x, y, z] = curveGeometry.data.position[0];
        return new THREE.Vector3(x, y, z);
      } else {
        console.warn(
          `Camera ${geometryKey}: curve key '${curveKey}' not found or invalid`
        );
        return fallback;
      }
    }
  };

  // Compute position from curve
  useEffect(() => {
    const point = computePointFromCurve(
      positionCurveKey,
      positionCurve,
      data.position_progress,
      new THREE.Vector3(0, 0, 10)
    );
    setComputedPosition(point);
  }, [positionCurve, positionCurveKey, data.position_progress, geometries, geometryKey]);

  // Compute target from curve
  useEffect(() => {
    const point = computePointFromCurve(
      targetCurveKey,
      targetCurve,
      data.target_progress,
      new THREE.Vector3(0, 0, 0)
    );
    setComputedTarget(point);
  }, [targetCurve, targetCurveKey, data.target_progress, geometries, geometryKey]);

  // Update scene camera if this camera is attached
  useEffect(() => {
    if (!isAttached) return;

    sceneCamera.position.copy(computedPosition);
    sceneCamera.up.set(data.up[0], data.up[1], data.up[2]);
    sceneCamera.lookAt(computedTarget);

    // Update projection properties
    let needsProjectionUpdate = false;

    if ('fov' in sceneCamera) {
      (sceneCamera as THREE.PerspectiveCamera).fov = data.fov;
      needsProjectionUpdate = true;
    }
    if ('near' in sceneCamera) {
      sceneCamera.near = data.near;
      needsProjectionUpdate = true;
    }
    if ('far' in sceneCamera) {
      sceneCamera.far = data.far;
      needsProjectionUpdate = true;
    }
    if ('zoom' in sceneCamera) {
      sceneCamera.zoom = data.zoom;
      needsProjectionUpdate = true;
    }

    // Update projection matrix if any projection property changed
    if (needsProjectionUpdate) {
      sceneCamera.updateProjectionMatrix();
    }
  }, [isAttached, computedPosition, computedTarget, data.up, data.fov, data.near, data.far, data.zoom, sceneCamera]);

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
              3
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
