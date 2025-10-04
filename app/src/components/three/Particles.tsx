import * as THREE from 'three';
import { useQueries } from '@tanstack/react-query';
import { getFrameDataOptions } from '../../hooks/useTrajectoryData';
import { useAppStore } from '../../store';
import { useRef, useMemo } from 'react';
import { useExtensionData } from '../../hooks/useSchemas';
import type { Representation } from '../../types/room-config';
import { useFrame } from '@react-three/fiber';

interface SphereProps {
  positionKey: string;
  colorKey: string;
  radiusKey: string;
}

// Reusable vectors and objects to avoid creating them in the loop
const positionVec = new THREE.Vector3();
const scaleVec = new THREE.Vector3();
const matrix = new THREE.Matrix4();
const color = new THREE.Color();

export default function Sphere({ positionKey, colorKey, radiusKey }: SphereProps) {
  const instancedMeshRef = useRef<THREE.InstancedMesh | null>(null);
  const { currentFrame, roomId, clientId, userId, selection } = useAppStore();
  const lastGoodFrameData = useRef<any>(null);

  const selectionSet = useMemo(() => {
    return selection ? new Set(selection) : null;
  }, [selection]);

  const { data: representationSettings } = useExtensionData(
    roomId || '',
    userId || '',
    'settings',
    'representation'
  ) as { data: Representation | undefined };

  const particleResolution = representationSettings?.particle_resolution ?? 8;
  const material = representationSettings?.material ?? 'MeshStandardMaterial';
  const particleScale = representationSettings?.particle_scale ?? 1.0;

  const requiredKeys = useMemo(() => [positionKey, colorKey, radiusKey], [
    positionKey,
    colorKey,
    radiusKey,
  ]);

  const queries = useMemo(() => {
    if (!roomId) {
      return [];
    }
    return requiredKeys.map(key => getFrameDataOptions(roomId, currentFrame, key));
  }, [currentFrame, roomId, requiredKeys]);

  const queryResults = useQueries({ queries });

  const { frameData, isFetching } = useMemo(() => {
    const isFetching = queryResults.some(result => result.isFetching || result.isPlaceholderData);
    const allDataPresent = queryResults.every(result => result.data);

    if (!allDataPresent || queryResults.length === 0) {
      return { isFetching, frameData: null };
    }

    const firstSuccessfulResult = queryResults.find(result => result.isSuccess);
    const combinedData: { [key: string]: any } = {};
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
  }, [queryResults, requiredKeys]);

  const dataToRender = frameData || lastGoodFrameData.current;

  useFrame(() => {
    if (!instancedMeshRef.current || !dataToRender || isFetching) {
      return;
    }

    const mesh = instancedMeshRef.current;
    
    const positions = dataToRender[positionKey];
    const colors = dataToRender[colorKey];
    const radii = dataToRender[radiusKey];
    const { count } = dataToRender;

    if (!positions || !colors || !radii || !count) {
      return;
    }

    for (let i = 0; i < count; i++) {
      const i3 = i * 3;
      positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
      scaleVec.set(radii[i] * particleScale, radii[i] * particleScale, radii[i] * particleScale);

      matrix.identity().setPosition(positionVec).scale(scaleVec);
      mesh.setMatrixAt(i, matrix);

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

  if (!clientId || !roomId || !dataToRender) {
    return null;
  }

  if (frameData) {
    lastGoodFrameData.current = frameData;
  }

  const renderMaterial = () => {
    const commonProps = { color: "white", side: THREE.FrontSide };
    switch (material) {
        case 'MeshBasicMaterial': return <meshBasicMaterial {...commonProps} />;
        case 'MeshPhysicalMaterial': return <meshPhysicalMaterial {...commonProps} roughness={0.3} reflectivity={0.4} />;
        case 'MeshToonMaterial': return <meshToonMaterial {...commonProps} />;
        case 'MeshStandardMaterial': default: return <meshStandardMaterial {...commonProps} />;
    }
  };

  return (
    <instancedMesh
      key={dataToRender.count} // Unique key based on count
      ref={instancedMeshRef}
      args={[undefined, undefined, dataToRender.count]}
    >
      <sphereGeometry args={[1, particleResolution, particleResolution]} />
      {renderMaterial()}
    </instancedMesh>
  );
}