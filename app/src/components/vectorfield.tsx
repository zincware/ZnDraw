import React from "react";
import { Canvas } from "@react-three/fiber";
import { Line } from "@react-three/drei";
import * as THREE from "three";

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
          <Line points={vector} color="blue" lineWidth={2} />
          {showArrows && vector.length > 1 && (
            <ArrowHelper
              start={vector[vector.length - 2]}
              end={vector[vector.length - 1]}
            />
          )}
        </React.Fragment>
      ))}
    </>
  );
};

interface ArrowHelperProps {
  start: [number, number, number];
  end: [number, number, number];
}

const ArrowHelper: React.FC<ArrowHelperProps> = ({ start, end }) => {
  const dir = new THREE.Vector3(
    end[0] - start[0],
    end[1] - start[1],
    end[2] - start[2],
  ).normalize();
  const length = new THREE.Vector3(
    end[0] - start[0],
    end[1] - start[1],
    end[2] - start[2],
  ).length();
  const arrowHelper = new THREE.ArrowHelper(
    dir,
    new THREE.Vector3(...start),
    length,
    0x0000ff,
  );

  return <primitive object={arrowHelper} />;
};

export default VectorField;
