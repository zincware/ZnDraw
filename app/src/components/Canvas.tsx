import * as THREE from 'three';
import { useQueries } from '@tanstack/react-query';
import { getFrameDataOptions } from '../hooks/useTrajectoryData';
import { useAppStore } from '../store';
import { useRef, useEffect, useMemo } from 'react';
import { Canvas, useFrame, useThree } from '@react-three/fiber';
import { OrbitControls, ContactShadows } from '@react-three/drei';
import { useExtensionData } from '../hooks/useSchemas';
import type { Representation } from '../types/room-config';
import { useColorScheme } from '@mui/material/styles';

// Reusable THREE objects
const matrix = new THREE.Matrix4();
const color = new THREE.Color();
const positionVec = new THREE.Vector3();
const scaleVec = new THREE.Vector3();

function Instances() {
  const instancedMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const { currentFrame, roomId, clientId, userId, selection } = useAppStore();
  const lastGoodFrameData = useRef<any>(null);
  const progress = useRef(0);
  
  // Memoize the selection Set for performance
  const selectionSet = useMemo(() => {
    console.log('Current selection:', selection);
    return selection ? new Set(selection) : null;
  }, [selection]);

  // Fetch representation settings
  const { data: representationSettings } = useExtensionData(
    roomId || '',
    userId || '',
    'settings',
    'representation'
  ) as { data: Representation | undefined };

  const particleResolution = representationSettings?.particle_resolution ?? 8;
  const material = representationSettings?.material ?? 'MeshStandardMaterial';
  const particleScale = representationSettings?.particle_scale ?? 1.0;

  const requiredKeys = ['positions', 'colors', 'radii'];

  // All hooks must be called unconditionally at the top level.
  // Memoize the query options. If roomId is missing, this will return an empty array.
  const queries = useMemo(() => {
    // FIX 1: Handle the case where roomId is not yet available inside the hook.
    if (!roomId) {
      return [];
    }
    return requiredKeys.map(key => getFrameDataOptions(roomId, currentFrame, key));
  }, [currentFrame, roomId]);

  // useQueries is safe to call with an empty array.
  const queryResults = useQueries({ queries });

  // This dependency array is better. It only changes when the data objects themselves change.
  const queryData = queryResults.map(q => q.data);
  const { frameData, isFetching } = useMemo(() => {
    const isFetching = queryResults.some(result => result.isFetching || result.isPlaceholderData);
    const allDataPresent = queryResults.every(result => result.data);

    // If queries haven't run or are still fetching, there's no data.
    if (!allDataPresent || queryResults.length === 0) {
      return { isFetching, frameData: null };
    }

    const firstSuccessfulResult = queryResults.find(result => result.isSuccess);
    const combinedData = {};
    let isComplete = true;

    for (let i = 0; i < requiredKeys.length; i++) {
      const key = requiredKeys[i];
      const result = queryResults[i];
      if (result?.data?.data) {
        combinedData[key] = result.data.data;
      } else {
        isComplete = false;
        break;
      }
    }

    if (!isComplete) {
      return { isFetching, frameData: null };
    }

    combinedData.count = firstSuccessfulResult?.data?.shape[0] || 0;
    return { isFetching, frameData: combinedData };
  }, [queryData]); // Depend on the array of data objects.

  // This hook can now run safely on every render.
  useEffect(() => {
    // Reset progress when new, valid data arrives.
    if (frameData) {
      progress.current = 0;
    }
  }, [frameData]);

  const dataToRender = frameData || lastGoodFrameData.current;

  // This hook also runs safely on every render.
  useFrame(() => {
    if (!instancedMeshRef.current || !dataToRender || isFetching) {
      return;
    }

    const mesh = instancedMeshRef.current;
    const { positions, colors, radii, count } = dataToRender;

    if (!positions || !colors || !radii || !count) {
      return;
    }

    if (progress.current >= count) {
      return;
    }

    // TODO: remove the chunking, it makes things only slower!
    const chunkSize = 1000000;
    const end = Math.min(progress.current + chunkSize, count);

    for (let i = progress.current; i < end; i++) {
      const i3 = i * 3;
      positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
      scaleVec.set(radii[i] * particleScale, radii[i] * particleScale, radii[i] * particleScale);

      matrix.identity().setPosition(positionVec).scale(scaleVec);
      mesh.setMatrixAt(i, matrix);

      // If this particle is in the selection, color it pink, otherwise use original color
      if (selectionSet && selectionSet.has(i)) {
        color.setRGB(1.0, 0.75, 0.8); // Pink color
      } else {
        color.setRGB(colors[i3], colors[i3 + 1], colors[i3 + 2]);
      }
      mesh.setColorAt(i, color);
    }

    mesh.instanceMatrix.needsUpdate = true;
    if (mesh.instanceColor) {
      mesh.instanceColor.needsUpdate = true;
    }

    progress.current = end;
  });

  // FIX 2: The guard is now placed *after* all hooks have been called.
  if (!clientId || !roomId || !dataToRender) {
    // You can return a loader here for the initial state
    return null;
  }

  if (frameData) {
    lastGoodFrameData.current = frameData;
  }

  // Render the appropriate material based on settings
  const renderMaterial = () => {
    const commonProps = { color: "white", side: THREE.FrontSide };

    switch (material) {
      case 'MeshBasicMaterial':
        return <meshBasicMaterial {...commonProps} />;
      case 'MeshPhysicalMaterial':
        return <meshPhysicalMaterial {...commonProps} roughness={0.3} reflectivity={0.4} />;
      case 'MeshToonMaterial':
        return <meshToonMaterial {...commonProps} />;
      case 'MeshStandardMaterial':
      default:
        return <meshStandardMaterial {...commonProps} />;
    }
  };

  return (
    <instancedMesh
      key={dataToRender.count}
      ref={instancedMeshRef}
      args={[undefined, undefined, dataToRender.count]}
    >
      <sphereGeometry args={[1, particleResolution, particleResolution]} />
      {renderMaterial()}
    </instancedMesh>
  );
}

function CameraAttachedLight({ intensity = 1.0 }) {
  const { camera } = useThree();
  const lightRef = useRef();

  useFrame(() => {
    if (lightRef.current) {
      // Copy the camera's world position
      lightRef.current.position.copy(camera.position);
      // Optional: You could also make the light target a point in front of the camera
    }
  });

  return <directionalLight ref={lightRef} intensity={intensity} />;
}

// ... MyScene component remains the same ...
function MyScene() {

  const { roomId, userId } = useAppStore();
  const { mode } = useColorScheme();
  const { data: studioLightingSettings } = useExtensionData(
    roomId || '',
    userId || '',
    'settings',
    'studio_lighting' // Assuming you have a studio lighting settings model
  );

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
        // Set the neutral background color
        style={{ background: backgroundColor }}
      >
        {/* 1. Main light attached to the camera for clarity */}
        <CameraAttachedLight intensity={keyLightIntensity} />

        {/* 2. Soft ambient light to lift the shadows */}
        <ambientLight intensity={fillLightIntensity} />

        {/* 3. Rim light from the back-top-right to create highlights and define shape */}
        <directionalLight
          position={[10, 20, -20]} // Behind and above
          intensity={rimLightIntensity}
        />

        {/* Your scene content */}
        <Instances />

        {showContactShadow &&
          <ContactShadows
            position={[0, -15, 0]}
            opacity={0.5}
            scale={80}
            blur={2}
            far={30}
            color="#000000" // Ensure shadows are neutral
          />
        }

        <OrbitControls enableDamping={false} />
      </Canvas>
    </div>
  );
}

export default MyScene;