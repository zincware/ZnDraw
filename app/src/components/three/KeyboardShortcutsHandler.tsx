import { useEffect } from "react";
import { useThree } from "@react-three/fiber";
import { useQueryClient } from "@tanstack/react-query";
import { useAppStore } from "../../store";
import * as THREE from "three";
import { TypedArray } from "../../myapi/client";
import { OrbitControls as OrbitControlsImpl } from "three-stdlib";
import { isPositionStatic } from "../../utils/geometryEditing";

// Keys for dynamic data references
const positionsKey = "arrays.positions";
const connectivityKey = "info.connectivity";

/**
 * Component that handles 3D-related keyboard shortcuts.
 * Must be rendered inside the Canvas to access Three.js context.
 */
export const KeyboardShortcutsHandler = () => {
	const controls = useThree((state) => state.controls) as OrbitControlsImpl;
	const queryClient = useQueryClient();

	const {
		currentFrame,
		selections,
		geometries,
		updateSelectionForGeometry,
		toggleInfoBoxes,
		mode,
		enterDrawingMode,
		exitDrawingMode,
		enterEditingMode,
		exitEditingMode,
		cycleTransformMode,
		setEditingSelectedAxis,
		editingSelectedAxis,
		roomId,
		saveFrameEdits,
	} = useAppStore();

	useEffect(() => {
		const handleKeyDown = (event: KeyboardEvent) => {
			// Check if focus is on an input element
			const target = event.target as HTMLElement;
			const isInputFocused =
				target.tagName === "INPUT" ||
				target.tagName === "TEXTAREA" ||
				target.isContentEditable;

			if (isInputFocused) return;

			// Handle Ctrl/Cmd + A for select all
			if ((event.ctrlKey || event.metaKey) && event.key === "a") {
				event.preventDefault();
				// Retrieve the positions array from the query cache, handling both {positionsKey: number[]} and number[] formats
				let positionsData = queryClient.getQueryData<
					TypedArray | { [positionsKey]: TypedArray }
				>(["frame", roomId, currentFrame, positionsKey]);
				let positions: TypedArray | undefined;
				if (
					positionsData &&
					typeof positionsData === "object" &&
					positionsKey in positionsData
				) {
					positions = (positionsData as { [positionsKey]: TypedArray })[
						positionsKey
					];
				} else if (positionsData) {
					positions = positionsData as TypedArray;
				}

				if (positions) {
					const count = positions.length / 3; // Assuming 3 components per vertex
					const allIndices = Array.from({ length: count }, (_, i) => i);
					updateSelectionForGeometry("particles", allIndices);
				}
				return;
			}
			// Handle c for center camera
			if (event.key === "c" || event.key === "C") {
				event.preventDefault();

				// Check if there are any selections across all geometries
				const hasAnySelection = Object.values(selections).some(
					(indices) => indices && indices.length > 0,
				);

				let centroid: THREE.Vector3 | null = null;

				if (hasAnySelection) {
					// Calculate centroid of ALL selected elements across all geometries
					let sumX = 0,
						sumY = 0,
						sumZ = 0,
						totalCount = 0;

					// Get dynamic positions from query cache (used by geometries with position="arrays.positions")
					let dynamicPositions: TypedArray | undefined;
					const positionsData = queryClient.getQueryData<
						TypedArray | { [positionsKey]: TypedArray }
					>(["frame", roomId, currentFrame, positionsKey]);
					if (
						positionsData &&
						typeof positionsData === "object" &&
						positionsKey in positionsData
					) {
						dynamicPositions = (
							positionsData as { [positionsKey]: TypedArray }
						)[positionsKey];
					} else if (positionsData) {
						dynamicPositions = positionsData as TypedArray;
					}

					// Get connectivity data for bonds (format: [atomA, atomB, bondOrder, ...])
					let connectivityData: TypedArray | undefined;
					const connData = queryClient.getQueryData<
						TypedArray | { [connectivityKey]: TypedArray }
					>(["frame", roomId, currentFrame, connectivityKey]);
					if (
						connData &&
						typeof connData === "object" &&
						connectivityKey in connData
					) {
						connectivityData = (connData as { [connectivityKey]: TypedArray })[
							connectivityKey
						];
					} else if (connData) {
						connectivityData = connData as TypedArray;
					}

					// Iterate over all geometry selections
					for (const [geometryKey, selectedIndices] of Object.entries(
						selections,
					)) {
						if (!selectedIndices || selectedIndices.length === 0) continue;

						const geometry = geometries[geometryKey];
						const positions = geometry?.data?.position;
						const connectivity = geometry?.data?.connectivity;

						// Check if this is a bond-like geometry (has connectivity)
						const hasDynamicConnectivity =
							typeof connectivity === "string" &&
							connectivity === connectivityKey &&
							connectivityData &&
							dynamicPositions;

						if (
							hasDynamicConnectivity &&
							connectivityData &&
							dynamicPositions
						) {
							// Bonds: add both atom positions for each selected bond
							for (const bondIdx of selectedIndices) {
								const connBaseIdx = bondIdx * 3; // [atomA, atomB, bondOrder]
								if (connBaseIdx + 1 < connectivityData.length) {
									const atomA = Number(connectivityData[connBaseIdx]);
									const atomB = Number(connectivityData[connBaseIdx + 1]);

									// Add atomA position
									const posA = atomA * 3;
									if (posA + 2 < dynamicPositions.length) {
										sumX += Number(dynamicPositions[posA]);
										sumY += Number(dynamicPositions[posA + 1]);
										sumZ += Number(dynamicPositions[posA + 2]);
										totalCount++;
									}

									// Add atomB position
									const posB = atomB * 3;
									if (posB + 2 < dynamicPositions.length) {
										sumX += Number(dynamicPositions[posB]);
										sumY += Number(dynamicPositions[posB + 1]);
										sumZ += Number(dynamicPositions[posB + 2]);
										totalCount++;
									}
								}
							}
						} else if (isPositionStatic(positions)) {
							// Static positions: number[][] format
							for (const idx of selectedIndices) {
								if (idx >= 0 && idx < positions.length) {
									const [x, y, z] = positions[idx];
									sumX += x;
									sumY += y;
									sumZ += z;
									totalCount++;
								}
							}
						} else if (
							typeof positions === "string" &&
							positions === positionsKey &&
							dynamicPositions
						) {
							// Dynamic positions: geometry references "arrays.positions" TypedArray
							for (const idx of selectedIndices) {
								const baseIdx = idx * 3;
								if (baseIdx + 2 < dynamicPositions.length) {
									sumX += Number(dynamicPositions[baseIdx]);
									sumY += Number(dynamicPositions[baseIdx + 1]);
									sumZ += Number(dynamicPositions[baseIdx + 2]);
									totalCount++;
								}
							}
						}
					}

					if (totalCount > 0) {
						centroid = new THREE.Vector3(
							sumX / totalCount,
							sumY / totalCount,
							sumZ / totalCount,
						);
					}
				}

				// Fallback: use all particles centroid if no selection
				if (!centroid) {
					let positionsData = queryClient.getQueryData<
						TypedArray | { [positionsKey]: TypedArray }
					>(["frame", roomId, currentFrame, positionsKey]);
					let positions: TypedArray | undefined;
					if (
						positionsData &&
						typeof positionsData === "object" &&
						positionsKey in positionsData
					) {
						positions = (positionsData as { [positionsKey]: TypedArray })[
							positionsKey
						];
					} else if (positionsData) {
						positions = positionsData as TypedArray;
					}

					if (positions && positions.length >= 3) {
						const count = positions.length / 3;
						let sumX = 0,
							sumY = 0,
							sumZ = 0;
						for (let i = 0; i < positions.length; i += 3) {
							sumX += Number(positions[i]);
							sumY += Number(positions[i + 1]);
							sumZ += Number(positions[i + 2]);
						}
						centroid = new THREE.Vector3(
							sumX / count,
							sumY / count,
							sumZ / count,
						);
					}
				}

				if (centroid && controls) {
					const camera = controls.object;
					const target = controls.target;

					// Move camera to look at centroid from a distance
					const direction = new THREE.Vector3();
					camera.getWorldDirection(direction);
					const distance = camera.position.distanceTo(target);
					const newPosition = centroid
						.clone()
						.add(direction.multiplyScalar(-distance));

					camera.position.copy(newPosition);
					controls.target.copy(centroid);
					controls.update();
				} else if (!controls) {
					console.warn("Camera controls not available");
				}
				return;
			}

			// Handle X/Y/Z for axis selection in editing mode
			if (
				mode === "editing" &&
				(event.key === "x" ||
					event.key === "X" ||
					event.key === "y" ||
					event.key === "Y" ||
					event.key === "z" ||
					event.key === "Z")
			) {
				event.preventDefault();
				const axis = event.key.toLowerCase() as "x" | "y" | "z";
				setEditingSelectedAxis(axis);
				return;
			}

			// Handle x/X for drawing mode (mutually exclusive with editing)
			if (event.key === "x" || event.key === "X") {
				event.preventDefault();
				if (mode === "drawing") {
					// Exit drawing mode
					exitDrawingMode();
				} else if (mode === "view") {
					// Enter drawing mode (only from view mode)
					enterDrawingMode(queryClient);
				}
				// If in editing mode, handled above
				return;
			}

			// Handle e/E for editing mode (mutually exclusive with drawing)
			if (event.key === "e" || event.key === "E") {
				event.preventDefault();
				if (mode === "editing") {
					// Exit editing mode
					exitEditingMode();
				} else if (mode === "view") {
					// Enter editing mode (only from view mode)
					enterEditingMode();
				}
				// If in drawing mode, do nothing (modes are exclusive)
				return;
			}

			// Handle t/T for cycling transform mode (only in editing mode)
			if ((event.key === "t" || event.key === "T") && mode === "editing") {
				event.preventDefault();
				void cycleTransformMode(); // async - fire and forget
				return;
			}

			// Handle s/S for saving frame edits (only in editing mode)
			if ((event.key === "s" || event.key === "S") && mode === "editing") {
				event.preventDefault();
				void saveFrameEdits(); // async - fire and forget
				return;
			}

			// Handle i/I for toggle info boxes
			if (event.key === "i" || event.key === "I") {
				event.preventDefault();
				toggleInfoBoxes();
				return;
			}

			// Handle r for rotate camera
			if (event.key === "r" || event.key === "R") {
				event.preventDefault();

				const clockwise = event.ctrlKey;

				if (!controls) {
					console.warn("Camera controls not available");
					return;
				}

				const camera = controls.object;
				const target = controls.target;

				// Get the current view direction
				const direction = new THREE.Vector3();
				camera.getWorldDirection(direction);

				// Calculate rotation axis (perpendicular to view direction)
				const axis = direction.clone().normalize();

				// Rotate camera position around target
				const position = camera.position.clone().sub(target);
				const angle =
					(event.shiftKey ? Math.PI / 8 : Math.PI / 32) * (clockwise ? -1 : 1);
				position.applyAxisAngle(axis, angle);
				camera.position.copy(position.add(target));

				// Update camera up vector
				const up = camera.up.clone();
				up.applyAxisAngle(axis, angle);
				camera.up.copy(up);

				controls.update();
				return;
			}
		};

		// Handle keyup to clear axis selection
		const handleKeyUp = (event: KeyboardEvent) => {
			if (
				mode === "editing" &&
				(event.key === "x" ||
					event.key === "X" ||
					event.key === "y" ||
					event.key === "Y" ||
					event.key === "z" ||
					event.key === "Z")
			) {
				const axis = event.key.toLowerCase() as "x" | "y" | "z";
				// Only clear if this is the axis we're tracking
				if (editingSelectedAxis === axis) {
					setEditingSelectedAxis(null);
				}
			}
		};

		window.addEventListener("keydown", handleKeyDown);
		window.addEventListener("keyup", handleKeyUp);
		return () => {
			window.removeEventListener("keydown", handleKeyDown);
			window.removeEventListener("keyup", handleKeyUp);
		};
	}, [
		controls,
		selections,
		geometries,
		updateSelectionForGeometry,
		queryClient,
		currentFrame,
		toggleInfoBoxes,
		mode,
		enterDrawingMode,
		exitDrawingMode,
		enterEditingMode,
		exitEditingMode,
		cycleTransformMode,
		setEditingSelectedAxis,
		editingSelectedAxis,
		roomId,
		saveFrameEdits,
	]);

	return null; // This component only handles keyboard events
};
