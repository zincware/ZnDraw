import React, { useMemo } from "react";
import * as THREE from "three";

interface ArrowProps {
  start: [number, number, number];
  end: [number, number, number];
  scale_vector_thickness?: boolean;
}

// TODO: fix start and end
// TODO: fix color
// TODO: find a good solution for scaling?
// TODO: provide an instanced arrow mesh version like: Arrows(start: list, end: list, colormap)
//   which does the color handling automatically and only use that one
const Arrow: React.FC<ArrowProps> = ({ start, end, scale_vector_thickness }) => {
  const cylinderRadius = 0.04;
  const cylinderHeight = 0.6;
  const coneRadius = 0.1;
  const coneHeight = 0.4;

  const { position, rotation, scale, color } = useMemo(() => {
    const startVector = new THREE.Vector3(...start);
    const endVector = new THREE.Vector3(...end);
    const direction = new THREE.Vector3().subVectors(endVector, startVector);
    const length = direction.length();
    let scale;
    if (scale_vector_thickness) {
      scale = new THREE.Vector3(length, length, length);
    } else {
      scale = new THREE.Vector3(1, length, 1);
    }

    const quaternion = new THREE.Quaternion().setFromUnitVectors(
      new THREE.Vector3(0, 1, 0),
      direction.clone().normalize(),
    );

    const eulerRotation = new THREE.Euler().setFromQuaternion(quaternion);
    const color = new THREE.Color().setHSL(length / 8 - 0.6, 1.0, 0.5);

    return {
      position: startVector,
      rotation: eulerRotation,
      scale: scale,
      color,
    };
  }, [start, end, scale_vector_thickness]);

  return (
    <group position={position} rotation={rotation} scale={scale}>
      <mesh>
        <cylinderGeometry
          args={[cylinderRadius, cylinderRadius, cylinderHeight]}
        />
        <meshStandardMaterial color={color} />
      </mesh>
      <mesh position={[0, (cylinderHeight + coneHeight) / 2, 0]}>
        <coneGeometry args={[coneRadius, coneHeight, 32]} />
        <meshStandardMaterial color={color} />
      </mesh>
    </group>
  );
};

export default Arrow;
