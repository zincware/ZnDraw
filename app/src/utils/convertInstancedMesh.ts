import * as THREE from "three";
import { BufferGeometryUtils } from "three/examples/jsm/Addons.js";

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
	material?: THREE.Material,
): THREE.Mesh {
	const baseGeometry = instancedMesh.geometry;
	const baseMaterial = instancedMesh.material;
	const count = instancedMesh.count;

	// Early return if no instances
	if (count === 0) {
		console.warn(
			"[convertInstancedMesh] Instance count is 0, returning empty mesh",
		);
		const emptyMaterial = Array.isArray(baseMaterial)
			? baseMaterial[0].clone()
			: baseMaterial.clone();
		return new THREE.Mesh(new THREE.BufferGeometry(), emptyMaterial);
	}

	// Validate base geometry
	if (!baseGeometry || !baseGeometry.attributes.position) {
		throw new Error("Base geometry is invalid or missing position attribute");
	}

	// Reusable temporary objects (avoids creating new objects in loop)
	const tempMatrix = new THREE.Matrix4();
	const tempColor = new THREE.Color();

	// Get base attributes
	const basePositions = baseGeometry.attributes.position;
	const baseNormals = baseGeometry.attributes.normal;
	const baseUvs = baseGeometry.attributes.uv;
	const baseIndices = baseGeometry.index;

	const baseVertexCount = basePositions.count;
	const totalVertexCount = baseVertexCount * count;
	const baseIndexCount = baseIndices ? baseIndices.count : 0;
	const totalIndexCount = baseIndexCount * count;

	// Create new attributes for merged geometry
	const mergedPositions = new Float32Array(totalVertexCount * 3);
	const mergedNormals = baseNormals
		? new Float32Array(totalVertexCount * 3)
		: null;
	const mergedUvs = baseUvs ? new Float32Array(totalVertexCount * 2) : null;
	const mergedColors = new Float32Array(totalVertexCount * 3);
	// Use Uint32Array if base geometry uses it OR if totalVertexCount exceeds Uint16 max (65535)
	const mergedIndices = baseIndices
		? new (baseIndices.array instanceof Uint32Array || totalVertexCount > 65535
				? Uint32Array
				: Uint16Array)(totalIndexCount)
		: null;

	const hasInstanceColor = !!instancedMesh.instanceColor;

	// Transform and merge
	for (let i = 0; i < count; i++) {
		// Get instance matrix
		instancedMesh.getMatrixAt(i, tempMatrix);

		// Get instance color
		if (hasInstanceColor) {
			instancedMesh.getColorAt(i, tempColor);
		} else {
			tempColor.set(1, 1, 1); // Default white
		}

		const vertexOffset = i * baseVertexCount;
		const indexOffset = i * baseIndexCount;

		// Process vertices
		for (let v = 0; v < baseVertexCount; v++) {
			const targetIdx = (vertexOffset + v) * 3;
			const srcIdx = v * 3;

			// Position transform
			const x = basePositions.getX(v);
			const y = basePositions.getY(v);
			const z = basePositions.getZ(v);

			// Manual matrix multiplication for performance (avoid Vector3 allocation)
			const e = tempMatrix.elements;
			const w = 1 / (e[3] * x + e[7] * y + e[11] * z + e[15]);

			mergedPositions[targetIdx] = (e[0] * x + e[4] * y + e[8] * z + e[12]) * w;
			mergedPositions[targetIdx + 1] =
				(e[1] * x + e[5] * y + e[9] * z + e[13]) * w;
			mergedPositions[targetIdx + 2] =
				(e[2] * x + e[6] * y + e[10] * z + e[14]) * w;

			// Normal transform (if available)
			if (mergedNormals && baseNormals) {
				const nx = baseNormals.getX(v);
				const ny = baseNormals.getY(v);
				const nz = baseNormals.getZ(v);

				// Transform normal by rotation part of matrix (upper 3x3)
				// Note: For non-uniform scale, we technically need inverse-transpose,
				// but for uniform scale (common in particles) this is sufficient.
				// For perfect correctness with non-uniform scale, we'd need normal matrix.
				// Given performance constraints, direct rotation is often acceptable approximation here.
				mergedNormals[targetIdx] = e[0] * nx + e[4] * ny + e[8] * nz;
				mergedNormals[targetIdx + 1] = e[1] * nx + e[5] * ny + e[9] * nz;
				mergedNormals[targetIdx + 2] = e[2] * nx + e[6] * ny + e[10] * nz;
			}

			// UV copy (if available)
			if (mergedUvs && baseUvs) {
				const targetUvIdx = (vertexOffset + v) * 2;
				mergedUvs[targetUvIdx] = baseUvs.getX(v);
				mergedUvs[targetUvIdx + 1] = baseUvs.getY(v);
			}

			// Color copy
			mergedColors[targetIdx] = tempColor.r;
			mergedColors[targetIdx + 1] = tempColor.g;
			mergedColors[targetIdx + 2] = tempColor.b;
		}

		// Copy indices (offset by vertex count)
		if (mergedIndices && baseIndices) {
			for (let j = 0; j < baseIndexCount; j++) {
				mergedIndices[indexOffset + j] = baseIndices.getX(j) + vertexOffset;
			}
		}
	}

	// Create merged buffer geometry
	const mergedGeometry = new THREE.BufferGeometry();
	mergedGeometry.setAttribute(
		"position",
		new THREE.BufferAttribute(mergedPositions, 3),
	);
	if (mergedNormals)
		mergedGeometry.setAttribute(
			"normal",
			new THREE.BufferAttribute(mergedNormals, 3),
		);
	if (mergedUvs)
		mergedGeometry.setAttribute("uv", new THREE.BufferAttribute(mergedUvs, 2));
	mergedGeometry.setAttribute(
		"color",
		new THREE.BufferAttribute(mergedColors, 3),
	);
	if (mergedIndices)
		mergedGeometry.setIndex(new THREE.BufferAttribute(mergedIndices, 1));

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
		const clonedMaterial = Array.isArray(baseMaterial)
			? baseMaterial[0].clone()
			: baseMaterial.clone();

		// Enable vertex colors on the material
		if ("vertexColors" in clonedMaterial) {
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
			mesh.material.forEach((m) => m.dispose());
		} else {
			mesh.material.dispose();
		}
	}
}
