import React, { useRef, useState } from 'react';
import { Canvas } from "@react-three/fiber";
import { Line } from "@react-three/drei";
import * as THREE from "three";

interface ArrowProps {
  start: [number, number, number];
  end: [number, number, number];
}

const Arrow: React.FC<ArrowProps> = ({ start, end}) => {
  const cylinderRadius = 0.07;
  const cylinderHeight = 0.6;
  const coneRadius = 0.14;
  const coneHeight = 0.4;

  const rotation = new THREE.Euler().setFromQuaternion(
    new THREE.Quaternion().setFromUnitVectors(
      new THREE.Vector3(0, 1, 0),
      new THREE.Vector3(end[0] - start[0], end[1] - start[1], end[2] - start[2]),
    ),
  );

  const scale = new THREE.Vector3(start[0] - end[0], start[1] - end[1], start[2] - end[2]).length();

  return (
    <group position={start} rotation={rotation} scale={scale}>
      <mesh>
        <cylinderGeometry args={[cylinderRadius, cylinderRadius, cylinderHeight]} />
        <meshStandardMaterial color={'hotpink'} />
      </mesh>
      <mesh position={[0, (cylinderHeight + coneHeight ) / 2, 0]}>
        <coneGeometry args={[coneRadius, coneHeight, 32]} />
        <meshStandardMaterial color={'hotpink'} />
      </mesh>
    </group>
  );
};


interface VectorFieldProps {
  vectors: [number, number, number][][];
  showArrows?: boolean;
}

const VectorField: React.FC<VectorFieldProps> = ({
  vectors,
  showArrows = true,
}) => {
  return (
    <>
      {vectors.map((vector, index) => (
        <React.Fragment key={index}>
          {/* <Line points={vector} color="blue" lineWidth={2} /> */}
          {showArrows && vector.length > 1 && (
            <Arrow
              start={vector[vector.length - 2]}
              end={vector[vector.length - 1]}
            />
          )}
        </React.Fragment>
      ))}
    </>
  );
};


export default VectorField;
