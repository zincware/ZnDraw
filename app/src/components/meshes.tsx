import React, { useMemo } from "react";
import * as THREE from "three";
import { BufferGeometryUtils } from "three/examples/jsm/Addons.js";
import { interpolateColor, HSLColor, ColorRange } from "./utils";

function createArrowMesh() {
  const cylinderRadius = 0.04;
  const cylinderHeight = 0.6;
  const coneRadius = 0.1;
  const coneHeight = 0.4;

  const cylinderGeometry = new THREE.CylinderGeometry(
    cylinderRadius,
    cylinderRadius,
    cylinderHeight,
    32,
  );
  const coneGeometry = new THREE.ConeGeometry(coneRadius, coneHeight, 32);

  cylinderGeometry.translate(0, cylinderHeight / 2, 0);
  coneGeometry.translate(0, cylinderHeight + coneHeight / 2, 0);

  const arrowGeometry = BufferGeometryUtils.mergeGeometries([
    cylinderGeometry,
    coneGeometry,
  ]);

  return arrowGeometry;
}

interface ArrowsProps {
  start: number[][];
  end: number[][];
  scale_vector_thickness?: boolean;
  colormap: HSLColor[];
  colorrange: ColorRange;
  opacity?: number;
}

const Arrows: React.FC<ArrowsProps> = ({
  start,
  end,
  scale_vector_thickness,
  colormap,
  colorrange,
  opacity = 1.0,
}) => {
  const geometry = useMemo(() => createArrowMesh(), []);
  return (
    <>
      {start.map((s, i) => {
        const startVector = new THREE.Vector3(...s);
        const endVector = new THREE.Vector3(...end[i]);
        const direction = new THREE.Vector3().subVectors(
          endVector,
          startVector,
        );
        const length = direction.length();

        const scale = scale_vector_thickness
          ? new THREE.Vector3(length, length, length)
          : new THREE.Vector3(1, length, 1);
        const quaternion = new THREE.Quaternion().setFromUnitVectors(
          new THREE.Vector3(0, 1, 0),
          direction.clone().normalize(),
        );
        const eulerRotation = new THREE.Euler().setFromQuaternion(quaternion);
        const color = interpolateColor(colormap, colorrange, length);

        return (
          <mesh
            geometry={geometry}
            position={startVector}
            rotation={eulerRotation}
            scale={scale}
          >
            <meshStandardMaterial color={color} transparent opacity={opacity} />
          </mesh>
        );
      })}
    </>
  );
};

export default Arrows;
