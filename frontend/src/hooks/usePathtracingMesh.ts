import { type RefObject, useCallback, useEffect, useRef } from "react";
import type * as THREE from "three";
import { useAppStore } from "../store";
import {
	convertInstancedMeshToMerged,
	disposeMesh,
} from "../utils/convertInstancedMesh";

/**
 * Hook that manages the merged mesh for GPU pathtracing.
 *
 * GPU pathtracer doesn't support instanced meshes, so we convert them to a
 * single merged mesh with vertex colors. This hook handles the conversion,
 * disposal, and update signaling.
 *
 * CRITICAL: This hook uses refs and manual scene graph management (not React state)
 * to ensure precise timing control. The merged mesh is added to the parent group
 * synchronously, guaranteeing it's in the scene before the pathtracer updates.
 *
 * @param parentGroupRef - Ref to the parent group where merged mesh will be added
 * @param mainMeshRef - Ref to the instanced mesh to convert
 * @param pathtracingEnabled - Whether pathtracing mode is active
 * @returns updateMergedMesh function to trigger conversion
 */
export function usePathtracingMesh(
	parentGroupRef: RefObject<THREE.Group | null>,
	mainMeshRef: RefObject<THREE.InstancedMesh | null>,
	pathtracingEnabled: boolean,
): (geometry: THREE.BufferGeometry) => void {
	const requestPathtracingUpdate = useAppStore(
		(state) => state.requestPathtracingUpdate,
	);

	// Ref to track the merged mesh - NO React state to avoid timing issues
	const mergedMeshRef = useRef<THREE.Mesh | null>(null);

	// Function to convert and update the merged mesh
	// This is synchronous - mesh is added to scene immediately
	const updateMergedMesh = useCallback(
		(geometry: THREE.BufferGeometry) => {
			if (
				!pathtracingEnabled ||
				!mainMeshRef.current ||
				!parentGroupRef.current
			) {
				return;
			}

			// Remove and dispose old mesh from scene
			if (mergedMeshRef.current) {
				parentGroupRef.current.remove(mergedMeshRef.current);
				disposeMesh(mergedMeshRef.current);
				mergedMeshRef.current = null;
			}

			// Convert instanced mesh to merged mesh with vertex colors
			const newMergedMesh = convertInstancedMeshToMerged(
				mainMeshRef.current,
				geometry,
			);

			// Add new mesh to scene SYNCHRONOUSLY
			parentGroupRef.current.add(newMergedMesh);
			mergedMeshRef.current = newMergedMesh;

			// Signal pathtracer to rebuild BVH
			// The mesh is already in the scene, so update() will see it
			requestPathtracingUpdate();
		},
		[pathtracingEnabled, mainMeshRef, parentGroupRef, requestPathtracingUpdate],
	);

	// Cleanup when pathtracing is disabled
	useEffect(() => {
		if (
			!pathtracingEnabled &&
			mergedMeshRef.current &&
			parentGroupRef.current
		) {
			parentGroupRef.current.remove(mergedMeshRef.current);
			disposeMesh(mergedMeshRef.current);
			mergedMeshRef.current = null;
		}
	}, [pathtracingEnabled, parentGroupRef]);

	// Cleanup on unmount
	useEffect(() => {
		return () => {
			if (mergedMeshRef.current) {
				// Try to remove from parent if still attached
				if (mergedMeshRef.current.parent) {
					mergedMeshRef.current.parent.remove(mergedMeshRef.current);
				}
				disposeMesh(mergedMeshRef.current);
				mergedMeshRef.current = null;
			}
		};
	}, []);

	return updateMergedMesh;
}
