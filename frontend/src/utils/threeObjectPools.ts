import * as THREE from "three";

/**
 * Reusable THREE.js objects to avoid creating new instances in render loops.
 * These objects reduce garbage collection pressure and improve performance.
 *
 * IMPORTANT: These objects are shared across all geometry components.
 * Always clone or copy values before storing them permanently.
 *
 * Usage example:
 * ```typescript
 * _vec3.set(x, y, z);
 * mesh.position.copy(_vec3); // Copy, don't assign directly!
 * ```
 */

// Vector3 objects for position, direction, scale calculations
export const _vec3 = new THREE.Vector3();
export const _vec3_2 = new THREE.Vector3();
export const _vec3_3 = new THREE.Vector3();
export const _vec3_4 = new THREE.Vector3();

// Quaternion for rotation calculations
export const _quat = new THREE.Quaternion();
export const _quat2 = new THREE.Quaternion();

// Euler for rotation calculations
export const _euler = new THREE.Euler();

// Matrix for instance transformations
export const _matrix = new THREE.Matrix4();
export const _matrix2 = new THREE.Matrix4();

// Color for instance color calculations
export const _color = new THREE.Color();
