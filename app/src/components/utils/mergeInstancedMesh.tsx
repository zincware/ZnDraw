import * as THREE from 'three';
import * as BufferGeometryUtils from 'three/addons/utils/BufferGeometryUtils.js';
/**
 * Merges an InstancedMesh into a single Mesh.
 * 
 * @param {THREE.InstancedMesh} instancedMesh - The instanced mesh to merge.
 * @returns {THREE.Mesh} - A single mesh with merged geometry.
 */
export function mergeInstancedMesh(instancedMesh: THREE.InstancedMesh): THREE.Mesh {
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
    instanceGeometry.applyMatrix4(new THREE.Matrix4().compose(tempPosition, tempQuaternion, tempScale));

    // Add transformed geometry to the array for merging
    mergedGeometries.push(instanceGeometry);
  }

  // Merge all transformed geometries into one
  const mergedGeometry = BufferGeometryUtils.mergeGeometries(mergedGeometries, true);

  // Return a single mesh
  return new THREE.Mesh(mergedGeometry, material);
}