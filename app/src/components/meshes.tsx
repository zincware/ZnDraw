import React, { useMemo } from "react";
import * as THREE from "three";
import { interpolateColor, HSLColor, ColorRange } from "./utils";

interface ArrowProps {
  start: [number, number, number];
  end: [number, number, number];
  scale_vector_thickness?: boolean;
  colormap: HSLColor[];
  colorrange: ColorRange;
}

// TODO: fix start and end
// TODO: fix color
// TODO: find a good solution for scaling?
// TODO: provide an instanced arrow mesh version like: Arrows(start: list, end: list, colormap)
//   which does the color handling automatically and only use that one
const Arrow: React.FC<ArrowProps> = ({
  start,
  end,
  scale_vector_thickness,
  colormap,
  colorrange,
}) => {
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
    const color = interpolateColor(colormap, colorrange, length);
    console.log(
      "using color range",
      colorrange,
      "for length",
      length,
      "resulting in color",
      color,
    );

    return {
      position: startVector,
      rotation: eulerRotation,
      scale: scale,
      color,
    };
  }, [start, end, scale_vector_thickness, colormap, colorrange]);

  return (
    <group position={position} rotation={rotation} scale={scale}>
      <mesh position={[0, cylinderHeight / 2, 0]}>
        <cylinderGeometry
          args={[cylinderRadius, cylinderRadius, cylinderHeight]}
        />
        <meshStandardMaterial color={color} />
      </mesh>
      <mesh position={[0, cylinderHeight + coneHeight / 2, 0]}>
        <coneGeometry args={[coneRadius, coneHeight, 32]} />
        <meshStandardMaterial color={color} />
      </mesh>
    </group>
  );
};

export default Arrow;
