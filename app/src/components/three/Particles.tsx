import * as THREE from 'three';
import { useQueries } from '@tanstack/react-query';
import { getFrameDataOptions } from '../../hooks/useTrajectoryData';
import { useAppStore } from '../../store';
import { useRef, useEffect, useMemo } from 'react';
import { useExtensionData } from '../../hooks/useSchemas';
import type { Representation } from '../../types/room-config';
import { useFrame } from '@react-three/fiber';

// Reusable vectors to avoid creating them in the loop
const positionVec = new THREE.Vector3();
const scaleVec = new THREE.Vector3();

// Reusable THREE objects to avoid creating them in the loop
const matrix = new THREE.Matrix4();
const color = new THREE.Color();
export default function Particles() {
  const instancedMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const { currentFrame, roomId, clientId, userId, selection } = useAppStore();
  const lastGoodFrameData = useRef<any>(null);
  
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

  const requiredKeys = ['arrays.positions', 'arrays.colors', 'arrays.radii'];

  const queries = useMemo(() => {
    if (!roomId) {
      return [];
    }
    return requiredKeys.map(key => getFrameDataOptions(roomId, currentFrame, key));
  }, [currentFrame, roomId]);

  const queryResults = useQueries({ queries });

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

  const dataToRender = frameData || lastGoodFrameData.current;

  // This hook also runs safely on every render.
  useFrame(() => {
    if (!instancedMeshRef.current || !dataToRender || isFetching) {
      return;
    }

    const mesh = instancedMeshRef.current;
    const { "arrays.positions": positions, "arrays.colors": colors, "arrays.radii": radii, count } = dataToRender;

    if (!positions || !colors || !radii || !count) {
      return;
    }

    for (let i = 0; i < count; i++) {
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
