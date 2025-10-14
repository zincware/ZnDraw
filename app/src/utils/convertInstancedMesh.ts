import * as THREE from 'three';
import { BufferGeometryUtils } from 'three/examples/jsm/Addons.js';

/**
 * Converts an InstancedMesh to a single merged Mesh with vertex colors.
 * Required for GPU path tracing which does not support instancing.
 *
 * This is FAR more efficient than creating N individual meshes:
 * - Memory: 1 geometry + 1 material instead of N geometries + N materials
 * - Draw calls: 1 draw call instead of N draw calls
 * - Instance colors are baked into vertex colors
 *
 * @param instancedMesh - The instanced mesh to convert
 * @param material - Optional material to use (must support vertexColors)
 * @returns A single THREE.Mesh with merged geometry and vertex colors
 */
export function convertInstancedMeshToMerged(
  instancedMesh: THREE.InstancedMesh,
  material?: THREE.Material
): THREE.Mesh {
  const baseGeometry = instancedMesh.geometry;
  const baseMaterial = instancedMesh.material;
  const count = instancedMesh.count;

  // Early return if no instances
  if (count === 0) {
    console.warn('[convertInstancedMesh] Instance count is 0, returning empty mesh');
    const emptyMaterial = Array.isArray(baseMaterial) ? baseMaterial[0].clone() : baseMaterial.clone();
    return new THREE.Mesh(new THREE.BufferGeometry(), emptyMaterial);
  }

  // Validate base geometry
  if (!baseGeometry || !baseGeometry.attributes.position) {
    throw new Error('Base geometry is invalid or missing position attribute');
  }

  const tempMatrix = new THREE.Matrix4();
  const tempColor = new THREE.Color();

  // Array to hold all transformed geometries
  const geometries: THREE.BufferGeometry[] = [];

  // Temporary objects for matrix decomposition
  const tempPosition = new THREE.Vector3();
  const tempQuaternion = new THREE.Quaternion();
  const tempScale = new THREE.Vector3();

  for (let i = 0; i < count; i++) {
    // Get the transformation matrix for the current instance
    instancedMesh.getMatrixAt(i, tempMatrix);

    // Decompose the matrix into position, rotation (as quaternion), and scale
    tempMatrix.decompose(tempPosition, tempQuaternion, tempScale);

    // Clone the original geometry and apply the instance's transformation
    const instanceGeometry = baseGeometry.clone();
    instanceGeometry.applyMatrix4(
      new THREE.Matrix4().compose(tempPosition, tempQuaternion, tempScale)
    );

    // Bake instance color into vertex colors
    if (instancedMesh.instanceColor) {
      instancedMesh.getColorAt(i, tempColor);

      const colors: number[] = [];
      const vertexCount = instanceGeometry.attributes.position.count;
      for (let j = 0; j < vertexCount; j++) {
        colors.push(tempColor.r, tempColor.g, tempColor.b);
      }
      instanceGeometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    }

    // Add transformed geometry to the array for merging
    geometries.push(instanceGeometry);
  }

  // Merge all geometries into a single BufferGeometry
  // CRITICAL: Pass `true` as second parameter to create geometry groups
  let mergedGeometry: THREE.BufferGeometry;
  try {
    mergedGeometry = BufferGeometryUtils.mergeGeometries(geometries, true);
  } catch (error) {
    console.error('[convertInstancedMesh] Error merging geometries:', error);
    throw error;
  }

  // Dispose individual geometries after merging to free memory
  geometries.forEach(geo => geo.dispose());

  // CRITICAL: Compute bounding box and sphere for pathtracer!
  // The pathtracer needs these for ray intersection tests
  mergedGeometry.computeBoundingBox();
  mergedGeometry.computeBoundingSphere();

  // Use provided material or clone from base material with vertex colors enabled
  let finalMaterial: THREE.Material;

  if (material) {
    finalMaterial = material;
  } else {
    // Clone the base material and enable vertex colors
    const clonedMaterial = Array.isArray(baseMaterial) ? baseMaterial[0].clone() : baseMaterial.clone();

    // Enable vertex colors on the material
    if ('vertexColors' in clonedMaterial) {
      (clonedMaterial as any).vertexColors = true;
    }

    finalMaterial = clonedMaterial;
  }


  // Create the single merged mesh
  const mergedMesh = new THREE.Mesh(mergedGeometry, finalMaterial);
  return mergedMesh;
}

/**
 * Disposes of a mesh and its geometry/material.
 * Essential for preventing memory leaks when switching modes.
 */
export function disposeMesh(mesh: THREE.Mesh): void {
  if (mesh.geometry) {
    mesh.geometry.dispose();
  }
  if (mesh.material) {
    if (Array.isArray(mesh.material)) {
      mesh.material.forEach(m => m.dispose());
    } else {
      mesh.material.dispose();
    }
  }
}
