import * as THREE from 'three';
import { useQueries } from '@tanstack/react-query';
import { getFrameDataOptions } from '../hooks/useTrajectoryData';
import { useAppStore } from '../store';
import { useRef, useEffect, useMemo } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';

// Reusable THREE objects
const matrix = new THREE.Matrix4();
const color = new THREE.Color();
const positionVec = new THREE.Vector3();
const scaleVec = new THREE.Vector3();

function Instances() {
  const instancedMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const { currentFrame, roomId, clientId } = useAppStore();
  const lastGoodFrameData = useRef<any>(null);
  const progress = useRef(0);
  
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

    const chunkSize = 10000;
    const end = Math.min(progress.current + chunkSize, count);

    for (let i = progress.current; i < end; i++) {
      const i3 = i * 3;
      positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
      scaleVec.set(radii[i], radii[i], radii[i]);
      
      matrix.identity().setPosition(positionVec).scale(scaleVec);
      mesh.setMatrixAt(i, matrix);

      color.setRGB(colors[i3], colors[i3 + 1], colors[i3 + 2]);
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

  return (
    <instancedMesh
      key={dataToRender.count}
      ref={instancedMeshRef}
      args={[undefined, undefined, dataToRender.count]}
    >
      <sphereGeometry args={[1, 8, 8]} />
      <meshPhysicalMaterial color={"white"} />
    </instancedMesh>
  );
}

// ... MyScene component remains the same ...
function MyScene() {
  return (
    <div style={{ width: '100%', height: 'calc(100vh - 64px)' }}>
      <Canvas style={{ background: '#d6d6d6ff' }} >
        <OrbitControls />
        <ambientLight intensity={1.5} />
        <pointLight position={[10, 10, 10]} intensity={0.5} />
        <Instances />
      </Canvas>
    </div>
  );
}

export default MyScene;