import { useRef, useCallback } from "react";
import { Plane } from "@react-three/drei";
import { useFrame, useThree } from "@react-three/fiber";
import * as THREE from "three";
import { useAppStore } from "../../store";

/**
 * VirtualCanvas - An invisible plane that tracks mouse position in 3D space
 * during drawing mode. It catches pointer events when not hovering over
 * drawable geometries and provides the 3D cursor position for drawing.
 */
export default function VirtualCanvas() {
  const { camera } = useThree();
  // Use individual selectors to prevent unnecessary re-renders
  const mode = useAppStore((state) => state.mode);
  const setDrawingPointerPosition = useAppStore((state) => state.setDrawingPointerPosition);
  const drawingIsValid = useAppStore((state) => state.drawingIsValid);
  const planeRef = useRef<THREE.Mesh>(null);

  const isDrawing = mode === 'drawing';
  
  const distance = 10; // Fixed distance from camera

  // Update plane size and position every frame
  useFrame(() => {
    if (!planeRef.current || !isDrawing) return;

    // Calculate plane size based on camera FOV to fill the view
    // Type guard for PerspectiveCamera
    if (!('fov' in camera) || !('aspect' in camera)) return;
    
    const vFOV = THREE.MathUtils.degToRad(camera.fov);
    const height = 2 * Math.tan(vFOV / 2) * distance;
    const width = height * camera.aspect;
    
    planeRef.current.scale.set(width, height, 1);

    // Position plane in front of camera
    const direction = new THREE.Vector3(0, 0, -1).applyQuaternion(camera.quaternion);
    const position = direction.multiplyScalar(distance).add(camera.position);
    planeRef.current.position.copy(position);
    
    // Make plane face the camera
    planeRef.current.lookAt(camera.position);
  });

  const handlePointerMove = (e: any) => {
    if (!isDrawing) return;
    if (drawingIsValid) return;
    setDrawingPointerPosition(e.point);
  };

  // Don't render if not in drawing mode
  if (!isDrawing) return null;

  return (
    <Plane
      ref={planeRef}
      args={[10, 10, 1, 1]}
      onPointerMove={handlePointerMove}
    >
      <meshBasicMaterial
        color="gray"
        opacity={0}
        transparent
        wireframe
        side={THREE.DoubleSide}
      />
    </Plane>
  );
}
