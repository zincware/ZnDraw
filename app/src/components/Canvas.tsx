import * as THREE from 'three';
import { useRef } from 'react';
import { Canvas, useFrame, useThree } from '@react-three/fiber';
import { OrbitControls, ContactShadows } from '@react-three/drei';
import { useAppStore } from '../store';
import { useExtensionData } from '../hooks/useSchemas';
import { useColorScheme } from '@mui/material/styles';

// Import our new, self-contained components
import { SimulationCell } from './three/SimulationCell';
import Particles from './three/Particles';

// A small helper component to keep a light attached to the camera
function CameraAttachedLight({ intensity = 1.0 }) {
  const { camera } = useThree();
  const lightRef = useRef<THREE.DirectionalLight>(null!);

  useFrame(() => {
    if (lightRef.current) {
      lightRef.current.position.copy(camera.position);
    }
  });

  return <directionalLight ref={lightRef} intensity={intensity} />;
}

// The main scene component
function MyScene() {
  const { roomId, userId } = useAppStore();
  const { mode } = useColorScheme();

  const { data: studioLightingSettings } = useExtensionData(
    roomId || '',
    userId || '',
    'settings',
    'studio_lighting'
  );

  // Use settings or provide sensible defaults
  const keyLightIntensity = studioLightingSettings?.key_light_intensity ?? 1.2;
  const fillLightIntensity = studioLightingSettings?.fill_light_intensity ?? 0.3;
  const rimLightIntensity = studioLightingSettings?.rim_light_intensity ?? 1.5;
  const backgroundColor = studioLightingSettings?.background_color ?? (mode === 'dark' ? '#333840' : '#f5f5f5');
  const showContactShadow = studioLightingSettings?.contact_shadow ?? true;

  return (
    <div style={{ width: '100%', height: 'calc(100vh - 64px)' }}>
      <Canvas
        shadows
        camera={{ position: [-10, 10, 30], fov: 50 }}
        gl={{ antialias: true, toneMapping: THREE.ACESFilmicToneMapping }}
        style={{ background: backgroundColor }}
      >
        <CameraAttachedLight intensity={keyLightIntensity} />
        <ambientLight intensity={fillLightIntensity} />
        <directionalLight position={[10, 20, -20]} intensity={rimLightIntensity} />

        {/* Render our clean, refactored components */}
        <Particles />
        <SimulationCell />

        {showContactShadow &&
          <ContactShadows
            position={[0, -15, 0]}
            opacity={0.5}
            scale={80}
            blur={2}
            far={30}
            color="#000000"
          />
        }

        <OrbitControls enableDamping={false} />
      </Canvas>
    </div>
  );
}

export default MyScene;
