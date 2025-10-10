import { useRef } from 'react';
import { useFrame, useThree } from "@react-three/fiber";
import * as THREE from 'three';

interface SceneLightingProps {
 ambient_light: number,
 key_light: number,
 fill_light: number,
 rim_light: number,
 hemisphere_light: number 
}

export default function SceneLighting(props: SceneLightingProps) {
  const rigRef = useRef<THREE.Group>(null!);
  const { camera } = useThree();

  useFrame(() => {
    if (rigRef.current) {
      rigRef.current.position.copy(camera.position);
      rigRef.current.quaternion.copy(camera.quaternion);
    }
  });

  return (
    <group ref={rigRef}>
      {/* Ambient light to fill shadows evenly */}
      <ambientLight intensity={props.ambient_light} color="#ffffff" />

      {/* Main "key" light: soft, frontal */}
      <directionalLight
        position={[5, 2, 8]}
        intensity={props.key_light}
      />

      {/* Fill light: adds soft illumination from below */}
      <directionalLight
        position={[-4, -1, 6]}
        intensity={props.fill_light}
        color="#a0c4ff"
      />

      {/* Rim/back light: subtle highlight separation */}
      <directionalLight
        position={[0, 0, -50]}
        intensity={props.rim_light}
        color="#fff0f5"
      />

      {/* Optional hemisphere light to simulate sky/ground reflection */}
      <hemisphereLight
        color="#e0f7ff"
        groundColor="#ffffff"
        intensity={props.hemisphere_light}
      />
    </group>
  );
}