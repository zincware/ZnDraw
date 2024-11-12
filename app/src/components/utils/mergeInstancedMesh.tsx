import * as THREE from "three";
import * as BufferGeometryUtils from "three/addons/utils/BufferGeometryUtils.js";
import { useMemo, useEffect, RefObject, useRef } from 'react';
import { useThree } from '@react-three/fiber';
import { usePathtracer } from "@react-three/gpu-pathtracer";

interface PathTracingSettings {
  enabled: boolean;
  [key: string]: any;
}

/**
 * Hook to create a merged mesh from an InstancedMesh, add it to the scene, and update it for path tracing.
 * 
 * @param {RefObject<THREE.InstancedMesh>} meshRef - React ref of the InstancedMesh.
 * @param {THREE.BufferGeometry} geometry - Geometry to use for the merged mesh.
 * @param {PathTracingSettings} settings - Settings for path tracing and merging.
 * @param {Array<unknown>} dependencies - Dependencies array for recalculating merged mesh.
 * 
 * @returns {THREE.Group | null} Merged mesh if created, otherwise null.
 */
export function useMergedMesh(
  meshRef: RefObject<THREE.InstancedMesh>,
  geometry: THREE.BufferGeometry,
  settings: PathTracingSettings,
  dependencies: unknown[] = []
): THREE.Group | null {
  const { scene } = useThree();
  const { update } = usePathtracer();

  const mergedMesh = useMemo(() => {
    if (meshRef.current?.instanceMatrix?.array?.length > 0 && settings?.enabled) {
      const singleMesh = splitInstancedMesh(meshRef.current, geometry, settings);
      return singleMesh;
    }
    return null;
  }, [geometry, settings, ...dependencies]);

  useEffect(() => {
    if (mergedMesh && settings?.enabled) {
      scene.add(mergedMesh);
      update();

      return () => {
        scene.remove(mergedMesh);
      };
    }
  }, [mergedMesh, settings]);

  return mergedMesh;
}

/**
 * Merges an InstancedMesh into a single Mesh.
 *
 * @param {THREE.InstancedMesh} instancedMesh - The instanced mesh to merge.
 * @returns {THREE.Mesh} - A single mesh with merged geometry.
 */
export function mergeInstancedMesh(
  instancedMesh: THREE.InstancedMesh,
): THREE.Mesh {
  const { count, geometry, material } = instancedMesh;
  const mergedGeometries = [];

  // Temporary objects to hold decomposed transformations
  const tempPosition = new THREE.Vector3();
  const tempQuaternion = new THREE.Quaternion();
  const tempScale = new THREE.Vector3();
  const tempMatrix = new THREE.Matrix4();

  for (let i = 0; i < count; i++) {
    // Get the transformation matrix for the current instance
    instancedMesh.getMatrixAt(i, tempMatrix);

    // Decompose the matrix into position, rotation (as quaternion), and scale
    tempMatrix.decompose(tempPosition, tempQuaternion, tempScale);

    // Clone the original geometry and apply the instance's transformation
    const instanceGeometry = geometry.clone();
    instanceGeometry.applyMatrix4(
      new THREE.Matrix4().compose(tempPosition, tempQuaternion, tempScale),
    );

    // Add transformed geometry to the array for merging
    mergedGeometries.push(instanceGeometry);
  }

  // Merge all transformed geometries into one
  const mergedGeometry = BufferGeometryUtils.mergeGeometries(
    mergedGeometries,
    true,
  );

  // Return a single mesh
  return new THREE.Mesh(mergedGeometry, material);
}

/**
 * Converts an InstancedMesh into a Group of individual meshes, each with its own color.
 *
 * @param {THREE.InstancedMesh} instancedMesh - The instanced mesh to split into individual meshes.
 * @returns {THREE.Group} - A group containing individual meshes.
 */
export function splitInstancedMesh(
  instancedMesh: THREE.InstancedMesh,
  geometry: THREE.BufferGeometry,
  pathTracingSettings: any,
): THREE.Group {
  const { count } = instancedMesh;
  const group = new THREE.Group();

  // Temporary objects to hold decomposed transformations
  const tempPosition = new THREE.Vector3();
  const tempQuaternion = new THREE.Quaternion();
  const tempScale = new THREE.Vector3();
  const tempMatrix = new THREE.Matrix4();

  // Accessing the instance color array if it exists
  const colorArray = instancedMesh.instanceColor?.array;

  for (let i = 0; i < count; i++) {
    // Get the transformation matrix for the current instance
    instancedMesh.getMatrixAt(i, tempMatrix);
    tempMatrix.decompose(tempPosition, tempQuaternion, tempScale);

    // Clone the geometry for each instance and apply the transformation
    const instanceGeometry = geometry.clone();
    // const instanceGeometry = new THREE.SphereGeometry(1, 32, 32);
    instanceGeometry.applyMatrix4(
      new THREE.Matrix4().compose(tempPosition, tempQuaternion, tempScale),
    );

    // Determine the color for this instance
    let color = new THREE.Color(0xffffff); // Default to white if no color attribute
    if (colorArray) {
      const r = colorArray[i * 3];
      const g = colorArray[i * 3 + 1];
      const b = colorArray[i * 3 + 2];
      color = new THREE.Color(r, g, b);
    }

    // Create a material with the instance's color
    const instanceMaterial = new THREE.MeshPhysicalMaterial({
      color: color,
      metalness: pathTracingSettings.metalness,
      roughness: pathTracingSettings.roughness,
      clearcoat: pathTracingSettings.clearcoat,
      clearcoatRoughness: pathTracingSettings.clearcoatRoughness,
    });

    // Create a mesh for this instance
    const instanceMesh = new THREE.Mesh(instanceGeometry, instanceMaterial);

    // Add the mesh to the group
    group.add(instanceMesh);
  }

  return group;
}
