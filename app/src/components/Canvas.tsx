import * as THREE from 'three';
import { useFrameData } from '../hooks/useTrajectoryData';
import { useAppStore } from '../store';
import { useRef, useEffect, useState, use } from 'react';
import { Canvas, useFrame } from '@react-three/fiber'; // Import Canvas and useFrame correctly
import { OrbitControls } from '@react-three/drei';

// Reusable objects to avoid creating new ones in the loop
const matrix = new THREE.Matrix4();
const color = new THREE.Color();

function Instances() {
  const instancedMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const { currentFrame } = useAppStore();


  
  const { data: frameData, isLoading } = useFrameData(currentFrame, ['positions', 'colors', 'radii']);

  // Ref to track the progress of the non-blocking loop
  const progress = useRef(0);
  
  // Reset progress when new frame data arrives
  useEffect(() => {
    progress.current = 0;
  }, [frameData]);

  useFrame(() => {
    if (!frameData || !instancedMeshRef.current) return;

    const mesh = instancedMeshRef.current;
    const { positions, colors, radii, count } = frameData;

    // If we've already processed all atoms for this frame, do nothing.
    if (progress.current >= count) {
      return;
    }

    // --- NON-BLOCKING LOOP ---
    // Process a chunk of atoms each frame instead of all at once.
    const chunkSize = 10000; // Tweak this number based on performance
    const end = Math.min(progress.current + chunkSize, count);
    
    for (let i = progress.current; i < end; i++) {
      const i3 = i * 3;
      
      // BUG FIX: Reset the matrix on each iteration
      matrix.identity();
      matrix.setPosition(positions[i3], positions[i3 + 1], positions[i3 + 2]);
      matrix.scale(new THREE.Vector3(radii[i], radii[i], radii[i]));
      mesh.setMatrixAt(i, matrix);
      
      color.setRGB(colors[i3], colors[i3 + 1], colors[i3 + 2]);
      mesh.setColorAt(i, color);
    }
    
    // Mark attributes for GPU update
    mesh.instanceMatrix.needsUpdate = true;
    if (mesh.instanceColor) {
      mesh.instanceColor.needsUpdate = true;
    }

    // Update our progress
    progress.current = end;
  });

  // TODO: keep showing the previous frame until the new one is ready instead of showing nothing
  if (isLoading) return null;
  if (!frameData) return null;

  return (
    <instancedMesh
      key={frameData.count} 
      ref={instancedMeshRef}
      args={[undefined, undefined, frameData.count]} 
    >
      <sphereGeometry args={[1, 8, 8]} />
      <meshPhysicalMaterial color={0xffffff} /> 
    </instancedMesh>
  );
}


function MyScene() {
  const { currentFrame } = useAppStore();

  return (
    <div style={{ width: '100%', height: 'calc(100vh - 64px)', overflow: 'hidden', userSelect: 'none', pointerEvents: 'none'}}>
      <Canvas style={{ background: '#d6d6d6ff', pointerEvents: 'auto' }} >
        <OrbitControls />
        <ambientLight intensity={1} />
        <Instances />
      </Canvas>
    </div>
  );
}

export default MyScene;