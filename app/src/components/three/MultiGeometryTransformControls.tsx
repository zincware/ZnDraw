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
		geometries,
		selections,
		notifyEditingChange,
		showSnackbar,
		updateSelectionForGeometry,
		setEditingCombinedCentroid,
	} = useAppStore();

	// Virtual object at centroid for transform controls
	const virtualObjectRef = useRef<THREE.Object3D>(null);

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

	// Don't render if not editing or no valid selections
	if (mode !== "editing" || !centroid || !hasValidSelections) {
		return null;
	}

	return (
		<>
			{/* Virtual object at centroid */}
			<object3D ref={virtualObjectRef} position={centroid} />

			{/* Transform controls attached to virtual object */}
			{virtualObjectRef.current && (
				<TransformControls
					object={virtualObjectRef.current}
					mode="translate"
					onChange={handleTransformChange}
					onMouseDown={handleDragStart}
					onMouseUp={handleDragEnd}
				/>
			)}
		</>
	);
}
