import React, { useEffect, useMemo, useRef } from "react";
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
  const meshRef = useRef<THREE.InstancedMesh>(null);
  const materialRef = useRef<THREE.MeshStandardMaterial>(null);

  useEffect(() => {
    if (!meshRef.current) return;
    const matrix = new THREE.Matrix4();
    const up = new THREE.Vector3(0, 1, 0);
    const startVector = new THREE.Vector3();
    const endVector = new THREE.Vector3();
    const direction = new THREE.Vector3();
    const quaternion = new THREE.Quaternion();

    for (let i = 0; i < start.length; i++) {
      // const startVector = new THREE.Vector3(...start[i]);
      startVector.fromArray(start[i]);
      endVector.fromArray(end[i]);
      direction.subVectors(endVector, startVector);
      const length = direction.length();
      const color = interpolateColor(colormap, colorrange, length);

      const scale = scale_vector_thickness
        ? new THREE.Vector3(length, length, length)
        : new THREE.Vector3(1, length, 1);

      quaternion.setFromUnitVectors(
        up,
        direction.clone().normalize(),
      );
      matrix.makeRotationFromQuaternion(quaternion);
      matrix.setPosition(startVector);
      matrix.scale(scale);

      meshRef.current.setColorAt(i, color);
      meshRef.current.setMatrixAt(i, matrix);
    }
    meshRef.current.instanceMatrix.needsUpdate = true;
  }, [start, end, scale_vector_thickness]);

  useEffect(() => {
    if (!materialRef.current) return;
    materialRef.current.needsUpdate = true; // TODO: check for particles as well
    if (!meshRef.current) return;
    if (!meshRef.current.instanceColor) return;
    meshRef.current.instanceColor.needsUpdate = true;
  }, [colormap, colorrange]);

  return (
    <instancedMesh ref={meshRef} args={[undefined, undefined, start.length]}>
      <bufferGeometry attach="geometry" {...geometry} />
      <meshStandardMaterial
        ref={materialRef}
        attach="material"
        transparent
        opacity={opacity}
      />
    </instancedMesh>
  );
};

export default Arrows;
