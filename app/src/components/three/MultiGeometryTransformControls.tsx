import { useRef, useEffect, useState, useCallback, useMemo } from "react";
import * as THREE from "three";
import { throttle } from "lodash";
import { TransformControls } from "@react-three/drei";
import { useAppStore } from "../../store";
import {
	validateStaticPositions,
	calculateCombinedCentroid,
} from "../../utils/geometryEditing";

/**
 * MultiGeometryTransformControls provides a single transform control
 * at the centroid of all selected instances across geometries.
 *
 * Responsibilities:
 * - Calculate combined centroid from selections
 * - Render transform controls at centroid
 * - Broadcast transform matrix changes via callback system
 * - Validate that selections have static positions
 *
 * Individual geometry components are responsible for:
 * - Subscribing to transform changes
 * - Applying transforms to their own data
 * - Persisting their own updates
 */
export default function MultiGeometryTransformControls() {
	const {
		mode,
		transformMode,
		editingSelectedAxis,
		geometries,
		selections,
		notifyEditingChange,
		showSnackbar,
		updateSelectionForGeometry,
		setEditingCombinedCentroid,
	} = useAppStore();

	// Virtual object at centroid for transform controls
	const virtualObjectRef = useRef<THREE.Object3D>(null);

	// Track when the virtual object has mounted (needed for TransformControls)
	const [virtualObjectMounted, setVirtualObjectMounted] = useState(false);

	// Current centroid position
	const [centroid, setCentroid] = useState<[number, number, number] | null>(
		null,
	);
	const [hasValidSelections, setHasValidSelections] = useState(false);

	// Track geometries that have been validated to prevent repeated deselection
	const validatedGeometriesRef = useRef<Set<string>>(new Set());

	// Track if we're currently dragging to prevent snap-back
	const isDraggingRef = useRef(false);

	// Track the last selection key to detect when selection changes (not just positions)
	const lastSelectionKeyRef = useRef<string | null>(null);

	// Generate a stable key from selections (ignores position changes)
	const selectionKey = useMemo(() => {
		const entries = Object.entries(selections)
			.filter(([_, indices]) => indices && indices.length > 0)
			.map(([key, indices]) => `${key}:${indices.sort().join(",")}`)
			.sort()
			.join("|");
		return entries;
	}, [selections]);

	// Validate selections ONCE when entering edit mode - deselect invalid ones
	useEffect(() => {
		if (mode !== "editing") {
			validatedGeometriesRef.current.clear();
			return;
		}

		// Validate that all selected positions are static
		const { invalidGeometries } = validateStaticPositions(
			geometries,
			selections,
		);

		// Only deselect if we haven't already warned about these geometries
		const newInvalidGeometries = invalidGeometries.filter(
			(key) => !validatedGeometriesRef.current.has(key),
		);

		if (newInvalidGeometries.length > 0) {
			showSnackbar(
				`Deselected instances with dynamic positions in: ${newInvalidGeometries.join(", ")}`,
				"warning",
			);

			// Mark these as validated
			newInvalidGeometries.forEach((key) =>
				validatedGeometriesRef.current.add(key),
			);

			// Deselect invalid geometries
			for (const geometryKey of newInvalidGeometries) {
				updateSelectionForGeometry(geometryKey, []);
			}
		}
		// Note: showSnackbar and updateSelectionForGeometry are stable Zustand actions, omit from deps
		// eslint-disable-next-line react-hooks/exhaustive-deps
	}, [mode, geometries, selections]);

	// Calculate centroid and position virtual object
	// IMPORTANT: Only reposition when entering edit mode or selection changes, NOT during drag
	useEffect(() => {
		if (mode !== "editing") {
			setCentroid(null);
			setHasValidSelections(false);
			setVirtualObjectMounted(false);
			setEditingCombinedCentroid(null);
			lastSelectionKeyRef.current = null;
			return;
		}

		// Skip repositioning if we're actively dragging (prevents snap-back)
		if (isDraggingRef.current) {
			return;
		}

		// Check if selection actually changed (not just positions)
		const selectionChanged = lastSelectionKeyRef.current !== selectionKey;
		const isFirstInit = lastSelectionKeyRef.current === null;

		// Only reposition on selection change or first initialization
		if (!selectionChanged && !isFirstInit) {
			return;
		}

		// Validate that all selected positions are static
		const { validSelections } = validateStaticPositions(geometries, selections);

		// Calculate combined centroid
		const newCentroid = calculateCombinedCentroid(geometries, validSelections);
		if (!newCentroid) {
			setCentroid(null);
			setHasValidSelections(false);
			setEditingCombinedCentroid(null);
			return;
		}

		setCentroid(newCentroid);
		setHasValidSelections(true);
		setEditingCombinedCentroid(newCentroid);
		lastSelectionKeyRef.current = selectionKey;

		// Position virtual object at centroid
		if (virtualObjectRef.current) {
			virtualObjectRef.current.position.set(
				newCentroid[0],
				newCentroid[1],
				newCentroid[2],
			);
			virtualObjectRef.current.rotation.set(0, 0, 0);
			virtualObjectRef.current.scale.set(1, 1, 1);
		}
	}, [mode, geometries, selections, selectionKey, setEditingCombinedCentroid]);

	// Throttle the notification to geometry components
	// This limits how often geometries update their Zustand state
	// Leading edge = true ensures first update is immediate
	// Trailing edge = true ensures final position is captured
	const throttledNotify = useMemo(
		() =>
			throttle(
				(matrix: THREE.Matrix4) => {
					notifyEditingChange(matrix);
				},
				50, // 50ms = max 20 updates/second
				{ leading: true, trailing: true },
			),
		[notifyEditingChange],
	);

	// Cleanup throttle on unmount
	useEffect(() => {
		return () => {
			throttledNotify.cancel();
		};
	}, [throttledNotify]);

	// Handle transform changes - broadcast to subscribed components
	const handleTransformChange = useCallback(() => {
		if (!virtualObjectRef.current) {
			return;
		}

		// Get current transform matrix
		const matrix = new THREE.Matrix4();
		matrix.compose(
			virtualObjectRef.current.position,
			virtualObjectRef.current.quaternion,
			virtualObjectRef.current.scale,
		);

		// Use throttled notification to limit update frequency
		throttledNotify(matrix);
	}, [throttledNotify]);

	// Handle drag start - prevent snap-back during drag
	const handleDragStart = useCallback(() => {
		isDraggingRef.current = true;
	}, []);

	// Handle drag end - allow repositioning again
	const handleDragEnd = useCallback(() => {
		isDraggingRef.current = false;
	}, []);

	// Callback ref to track when object3D mounts/unmounts
	const handleVirtualObjectRef = useCallback(
		(object: THREE.Object3D | null) => {
			virtualObjectRef.current = object;
			setVirtualObjectMounted(!!object);
		},
		[],
	);

	// Handle arrow key transforms when axis is selected
	useEffect(() => {
		if (mode !== "editing" || !editingSelectedAxis || !virtualObjectRef.current)
			return;

		const handleKeyDown = (event: KeyboardEvent) => {
			if (event.key !== "ArrowUp" && event.key !== "ArrowDown") return;
			if (!virtualObjectRef.current || !editingSelectedAxis) return;

			event.preventDefault();

			const direction = event.key === "ArrowUp" ? 1 : -1;
			const obj = virtualObjectRef.current;

			// Step sizes (shift for larger steps)
			const translateStep = event.shiftKey ? 1.0 : 0.1;
			const rotateStep = event.shiftKey ? Math.PI / 4 : Math.PI / 36; // 45° or 5°
			const scaleFactor = event.shiftKey ? 1.5 : 1.1;

			const axisIndex = { x: 0, y: 1, z: 2 }[editingSelectedAxis];

			switch (transformMode) {
				case "translate": {
					const delta = [0, 0, 0];
					delta[axisIndex] = direction * translateStep;
					obj.position.x += delta[0];
					obj.position.y += delta[1];
					obj.position.z += delta[2];
					break;
				}
				case "rotate": {
					const euler = new THREE.Euler();
					euler.copy(obj.rotation);
					const rotations = [euler.x, euler.y, euler.z];
					rotations[axisIndex] += direction * rotateStep;
					obj.rotation.set(rotations[0], rotations[1], rotations[2]);
					break;
				}
				case "scale": {
					const factor = direction > 0 ? scaleFactor : 1 / scaleFactor;
					const scales = [obj.scale.x, obj.scale.y, obj.scale.z];
					scales[axisIndex] *= factor;
					obj.scale.set(scales[0], scales[1], scales[2]);
					break;
				}
			}

			// Broadcast the change via the transform notification system
			const matrix = new THREE.Matrix4();
			matrix.compose(obj.position, obj.quaternion, obj.scale);
			notifyEditingChange(matrix);
		};

		window.addEventListener("keydown", handleKeyDown);
		return () => window.removeEventListener("keydown", handleKeyDown);
	}, [mode, editingSelectedAxis, transformMode, notifyEditingChange]);

	// Don't render if not editing or no valid selections
	if (mode !== "editing" || !centroid || !hasValidSelections) {
		return null;
	}

	return (
		<>
			{/* Virtual object at centroid */}
			<object3D ref={handleVirtualObjectRef} position={centroid} />

			{/* Transform controls attached to virtual object */}
			{virtualObjectMounted && virtualObjectRef.current && (
				<TransformControls
					object={virtualObjectRef.current}
					mode={transformMode}
					onChange={handleTransformChange}
					onMouseDown={handleDragStart}
					onMouseUp={handleDragEnd}
				/>
			)}
		</>
	);
}
