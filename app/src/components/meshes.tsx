import React from "react";
import * as THREE from "three";

interface ArrowProps {
    start: [number, number, number];
    end: [number, number, number];
  }
  
  const Arrow: React.FC<ArrowProps> = ({ start, end }) => {
    const cylinderRadius = 0.07;
    const cylinderHeight = 0.6;
    const coneRadius = 0.14;
    const coneHeight = 0.4;
  
    const rotation = new THREE.Euler().setFromQuaternion(
      new THREE.Quaternion().setFromUnitVectors(
        new THREE.Vector3(0, 1, 0),
        new THREE.Vector3(
          end[0] - start[0],
          end[1] - start[1],
          end[2] - start[2],
        ),
      ),
    );
  
    const scale = new THREE.Vector3(
      start[0] - end[0],
      start[1] - end[1],
      start[2] - end[2],
    ).length();
    const color = new THREE.Color();
    color.setHSL(scale / 8 - 0.6, 1.0, 0.5);
  
    return (
      <group position={start} rotation={rotation} scale={scale}>
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