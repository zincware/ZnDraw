// CameraManager.js (or place it in the same file as MyScene)

import { useEffect } from 'react';
import { useThree } from '@react-three/fiber';

function CameraManager({ settings }) {
  const { camera } = useThree();

  useEffect(() => {
    if (settings) {
      // Update properties that can be changed on the fly
      camera.near = settings.near_plane;
      camera.far = settings.far_plane;

      // You could also update other properties like fov for a perspective camera
      // if (camera.isPerspectiveCamera) {
      //   camera.fov = settings.fov;
      // }
      
      // CRITICAL: This tells Three.js to re-calculate the camera's projection matrix
      // with the new values. Without this, you won't see any change.
      camera.updateProjectionMatrix();
    }
  }, [settings, camera]); // Re-run this effect when settings or the camera object changes

  return null; // This component doesn't render anything itself
}

export default CameraManager;