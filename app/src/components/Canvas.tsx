import * as THREE from 'three';
import { useQueries } from '@tanstack/react-query';
import { getFrameDataOptions } from '../hooks/useTrajectoryData';
import { useAppStore } from '../store';
import { useRef, useEffect, useMemo } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';

// Reusable THREE objects to avoid creating new ones in the render loop for performance
const matrix = new THREE.Matrix4();
const color = new THREE.Color();
const positionVec = new THREE.Vector3();
const scaleVec = new THREE.Vector3();

function Instances() {
  const instancedMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const { currentFrame, roomId } = useAppStore();
  const lastGoodFrameData = useRef<any>(null);
  const progress = useRef(0);

  // Guard against running hooks before roomId is available
  if (!roomId) {
    return null;
  }

  const requiredKeys = ['positions', 'colors', 'radii'];

  // Memoize the query options array to prevent TanStack Query from refetching continuously
  const queries = useMemo(() => {
    return requiredKeys.map(key => getFrameDataOptions(roomId, currentFrame, key));
  }, [currentFrame, roomId]);
  
  // https://tanstack.com/query/v5/docs/framework/react/guides/parallel-queries
  const queryResults = useQueries({ queries });

  // TODO: how to access the isFetching state in other components?
  // Either refactor or set in Zustand store
  const { frameData, isFetching } = useMemo(() => {
    // check for q.isPlaceholderData to see, if the data is still loading
    const isFetching = queryResults.some(result => result.isFetching || result.isPlaceholderData);
    const firstSuccessfulResult = queryResults.find(result => result.isSuccess);

    const allDataPresent = queryResults.every(result => result.data);

    if (!allDataPresent) {
      return { isFetching, frameData: null };
    }

    const combinedData = {};
    let isComplete = true; // Flag to ensure all data for a complete frame is present

    // Rely on the guaranteed order of useQueries results instead of searching
    for (let i = 0; i < requiredKeys.length; i++) {
      const key = requiredKeys[i];
      const result = queryResults[i];

      if (result?.data?.data) {
        combinedData[key] = result.data.data;
      } else {
        isComplete = false; // If any key is missing, the frame is not ready
        break;
      }
    }
    
    // Only return a valid frameData object if all parts were successfully combined
    if (!isComplete) {
        return { isFetching, frameData: null };
    }

    combinedData.count = firstSuccessfulResult.data?.shape[0] || 0;

    return { isFetching, frameData: combinedData };
  }, [queryResults.map(q => q.data)]);


  if (frameData) {
    lastGoodFrameData.current = frameData;
  }
  const dataToRender = frameData || lastGoodFrameData.current;

  // Reset the animation progress when a new frame's data is ready
  useEffect(() => {
    progress.current = 0;
  }, [dataToRender]);

  useFrame(() => {
    // Top-level guards for the render loop
    if (!frameData || !instancedMeshRef.current || isFetching) {
      return;
    }

    const mesh = instancedMeshRef.current;
    const { positions, colors, radii, count } = frameData;

    // Stronger guard to prevent crashes if a data array is missing
    if (!positions || !colors || !radii) {
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
      
      matrix.identity();
      matrix.setPosition(positionVec);
      matrix.scale(scaleVec);
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

  if (!dataToRender) {
    return null; // Or a full-screen loading spinner for the very first frame
  }

  return (
    <instancedMesh
      key={dataToRender.count} // Re-creates the mesh if atom count changes
      ref={instancedMeshRef}
      args={[undefined, undefined, dataToRender.count]}
    >
      <sphereGeometry args={[1, 8, 8]} />
      <meshPhysicalMaterial color={"white"} />
    </instancedMesh>
  );
}

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