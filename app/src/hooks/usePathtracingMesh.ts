import {
	useRef,
	useEffect,
	useState,
	useCallback,
	type RefObject,
} from "react";
import * as THREE from "three";
import { useAppStore } from "../store";
import {
	convertInstancedMeshToMerged,
	disposeMesh,
} from "../utils/convertInstancedMesh";

interface UsePathtracingMeshResult {
	/** The merged mesh for rendering, or null when disabled/empty */
	mergedMesh: THREE.Mesh | null;
	/** Call this after updating the instanced mesh to trigger conversion */
	updateMergedMesh: () => void;
}

/**
 * Hook that manages the merged mesh for GPU pathtracing.
 *
 * GPU pathtracer doesn't support instanced meshes, so we convert them to a
 * single merged mesh with vertex colors. This hook handles the conversion,
 * disposal, and update signaling.
 *
 * Usage:
 * ```tsx
 * const { mergedMesh, updateMergedMesh } = usePathtracingMesh(
 *   mainMeshRef,
 *   pathtracingEnabled
 * );
 *
 * // In your main data effect, at the end:
 * if (pathtracingEnabled && instanceCount > 0) {
 *   updateMergedMesh();
 * }
 *
 * // In JSX:
 * {pathtracingEnabled && mergedMesh && <primitive object={mergedMesh} />}
 * ```
 *
 * @param mainMeshRef - Ref to the instanced mesh to convert
 * @param pathtracingEnabled - Whether pathtracing mode is active
 * @returns Object with mergedMesh and updateMergedMesh function
 */
export function usePathtracingMesh(
	mainMeshRef: RefObject<THREE.InstancedMesh | null>,
	pathtracingEnabled: boolean,
): UsePathtracingMeshResult {
	const requestPathtracingUpdate = useAppStore(
		(state) => state.requestPathtracingUpdate,
	);

	// State to hold the merged mesh and trigger re-renders
	const [mergedMesh, setMergedMesh] = useState<THREE.Mesh | null>(null);
	// Ref to track for disposal without causing re-render loops
	const mergedMeshRef = useRef<THREE.Mesh | null>(null);

	// Function to convert and update the merged mesh
	const updateMergedMesh = useCallback(() => {
		if (!pathtracingEnabled || !mainMeshRef.current) {
			return;
		}

		// Dispose old mesh before creating new one
		if (mergedMeshRef.current) {
			disposeMesh(mergedMeshRef.current);
		}

		// Convert instanced mesh to merged mesh with vertex colors
		const newMergedMesh = convertInstancedMeshToMerged(mainMeshRef.current);
		mergedMeshRef.current = newMergedMesh;
		setMergedMesh(newMergedMesh);

		// Signal pathtracer to rebuild BVH
		requestPathtracingUpdate();
	}, [pathtracingEnabled, mainMeshRef, requestPathtracingUpdate]);

	// Cleanup when pathtracing is disabled
	useEffect(() => {
		if (!pathtracingEnabled && mergedMeshRef.current) {
			disposeMesh(mergedMeshRef.current);
			mergedMeshRef.current = null;
			setMergedMesh(null);
		}
	}, [pathtracingEnabled]);

	// Cleanup on unmount
	useEffect(() => {
		return () => {
			if (mergedMeshRef.current) {
				disposeMesh(mergedMeshRef.current);
				mergedMeshRef.current = null;
			}
		};
	}, []);

	return { mergedMesh, updateMergedMesh };
}
