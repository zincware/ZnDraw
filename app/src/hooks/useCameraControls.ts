import { useEffect, useState } from "react";

interface ControlsState {
  enabled: boolean;
  enablePan: boolean;
  enableRotate: boolean;
  enableZoom: boolean;
}

/**
 * Hook to determine camera control states based on camera attachments.
 *
 * Control states:
 * - No attachment: full controls
 * - Position curve only: rotate only (can look around from fixed position)
 * - Target curve only: pan and zoom (can orbit around fixed target)
 * - Both curves: no controls (fully constrained cinematic mode)
 *
 * @param attachedCameraKey - The key of the currently attached camera geometry
 * @param geometries - All geometries in the scene
 * @returns Control state configuration for OrbitControls
 */
export function useCameraControls(
  attachedCameraKey: string | null,
  geometries: Record<string, any>
): ControlsState {
  const [controlsState, setControlsState] = useState<ControlsState>({
    enabled: true,
    enablePan: true,
    enableRotate: true,
    enableZoom: true,
  });

  useEffect(() => {
    if (!attachedCameraKey) {
      // State 1: Not attached - full controls
      setControlsState({
        enabled: true,
        enablePan: true,
        enableRotate: true,
        enableZoom: true,
      });
      return;
    }

    const camera = geometries[attachedCameraKey];
    if (!camera || camera.type !== "Camera") {
      // Invalid attachment, fallback to full controls
      setControlsState({
        enabled: true,
        enablePan: true,
        enableRotate: true,
        enableZoom: true,
      });
      return;
    }

    // Check if curves are attached (curves are always used in new model)
    const positionConstrained = !!camera.data.position_curve_key;
    const targetConstrained = !!camera.data.target_curve_key;

    if (positionConstrained && targetConstrained) {
      // State 5: Fully constrained - disable all controls (cinematic mode)
      setControlsState({
        enabled: false,
        enablePan: false,
        enableRotate: false,
        enableZoom: false,
      });
    } else if (positionConstrained && !targetConstrained) {
      // State 3: Position locked - can rotate only (security camera on rail)
      setControlsState({
        enabled: true,
        enablePan: false,
        enableRotate: true,
        enableZoom: false,
      });
    } else if (!positionConstrained && targetConstrained) {
      // State 4: Target locked - can pan and zoom (follow moving object)
      setControlsState({
        enabled: true,
        enablePan: true,
        enableRotate: false,
        enableZoom: true,
      });
    } else {
      // State 2: Neither constrained - full controls
      setControlsState({
        enabled: true,
        enablePan: true,
        enableRotate: true,
        enableZoom: true,
      });
    }
  }, [attachedCameraKey, geometries]);

  return controlsState;
}
