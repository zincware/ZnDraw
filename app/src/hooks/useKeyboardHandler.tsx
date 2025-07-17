import { useEffect, useRef } from "react";
import { getCentroid } from "../components/particlesEditor";
import { useAppContext } from "../contexts/AppContext";
import * as THREE from "three";

export const useKeyboardHandler = () => {
	const {
		playing,
		setPlaying,
		step,
		setStep,
		length,
		selectedFrames,
		bookmarks,
		setBookmarks,
		currentFrame,
		selectedIds,
		points,
		setPoints,
		selectedPoint,
		setSelectedPoint,
		showParticleInfo,
		setShowParticleInfo,
		isDrawing,
		setIsDrawing,
		cameraAndControls,
		setCameraAndControls,
	} = useAppContext();

	// Ref to track last removal time for throttling
	const lastRemovalTimeRef = useRef<number>(0);
	const REMOVAL_THROTTLE_MS = 100; // Minimum time between removals

	useEffect(() => {
		const handleKeyDown = (event: KeyboardEvent) => {
			const canvasContainer = document.querySelector(".canvas-container");
			const isCanvasFocused =
				canvasContainer && canvasContainer.contains(document.activeElement);
			const isBodyFocused = document.activeElement === document.body;

			if (!isCanvasFocused && !isBodyFocused) {
				return;
			}

			if (event.key === " ") {
				event.preventDefault();
				event.stopPropagation();
				setPlaying((prevPlaying) => !prevPlaying);
				return;
			}

			if (event.key === "ArrowRight") {
				event.preventDefault();
				event.stopPropagation();
				setPlaying(false);
				if (event.shiftKey) {
					// Jump to the bookmark after the current step
					const bookmarkKeys = Object.keys(bookmarks)
						.map(Number)
						.sort((a, b) => a - b);
					const nextBookmark = bookmarkKeys.find((key) => key > step);

					if (nextBookmark !== undefined) {
						setStep(nextBookmark);
					}
				} else {
					if (selectedFrames.indices.size > 0 && selectedFrames.active) {
						const nextFrame = Array.from(selectedFrames.indices).find(
							(frame) => frame > step,
						);
						if (nextFrame) {
							setStep(nextFrame);
						} else {
							setStep(Math.min(...selectedFrames.indices));
						}
					} else {
						setStep((prevStep) => (prevStep + 1 < length ? prevStep + 1 : 0));
					}
				}
			} else if (event.key === "ArrowLeft") {
				event.preventDefault();
				event.stopPropagation();
				setPlaying(false);
				if (event.shiftKey) {
					// Jump to the bookmark before the current step
					const bookmarkKeys = Object.keys(bookmarks)
						.map(Number)
						.sort((a, b) => b - a);
					const previousBookmark = bookmarkKeys.find((key) => key < step);

					if (previousBookmark !== undefined) {
						setStep(previousBookmark);
					}
				} else {
					// Move to the previous step, or wrap around to the end
					// check if selectedFrames length is greater than 0, then only jump
					// between selectedFrames
					if (selectedFrames.indices.size > 0 && selectedFrames.active) {
						const previousFrame = Array.from(selectedFrames.indices)
							.reverse()
							.find((frame) => frame < step);
						if (previousFrame) {
							setStep(previousFrame);
						} else {
							setStep(Math.max(...selectedFrames.indices));
						}
					} else {
						setStep((prevStep) =>
							prevStep - 1 >= 0 ? prevStep - 1 : length - 1,
						);
					}
				}
			} else if (event.key === "b") {
				// Add current step to bookmarks
				const name = prompt("Enter bookmark name:");
				if (name) {
					setBookmarks({
						...bookmarks,
						[step]: name,
					});
				}
			} else if (event.key === "Delete" || event.key === "Backspace") {
				if (event.shiftKey) {
					// Shift + Backspace: Remove last point from line with throttling
					if (points.length > 0) {
						const now = Date.now();
						if (now - lastRemovalTimeRef.current >= REMOVAL_THROTTLE_MS) {
							event.preventDefault();
							event.stopPropagation();
							lastRemovalTimeRef.current = now;

							const updatedPoints = points.slice(0, -1);
							setPoints(updatedPoints);
							// If we removed the selected point, clear selection
							if (selectedPoint && points.length > 0) {
								const lastPoint = points[points.length - 1];
								if (selectedPoint.equals(lastPoint)) {
									setSelectedPoint(null);
								}
							}
						}
					}
				} else {
					// Regular Delete/Backspace: Delete selected atoms or selected point
					if (selectedIds.size > 0) {
						// Logic for deleting atoms would go here
						// TODO: Implement atom deletion logic
					}
					// Delete selected point
					if (selectedPoint) {
						const updatedPoints = points.filter(
							(point) => !point.equals(selectedPoint),
						);
						setPoints(updatedPoints);
						setSelectedPoint(null);
					}
				}
			} else if (event.key === "Escape") {
				// Clear selection
				setSelectedPoint(null);
			} else if (event.key === "c" || event.key === "C") {
				// Center camera on centroid of selected particles, or all particles if none selected
				const centroid =
					selectedIds.size > 0
						? getCentroid(currentFrame.positions, selectedIds)
						: getCentroid(currentFrame.positions, new Set());

				setCameraAndControls((prev: any) => ({
					...prev,
					target: centroid,
				}));
			} else if (event.key === "i" || event.key === "I") {
				// Toggle particle info display
				setShowParticleInfo(!showParticleInfo);
			} else if (event.key === "x" || event.key === "X") {
				// Toggle drawing mode
				setIsDrawing(!isDrawing);
			} else if (event.key === "o" || event.key === "O") {
				// Reset camera to original position
				event.preventDefault();
				event.stopPropagation();

				if (currentFrame.positions.length === 0) {
					return;
				}

				// Calculate the camera positions
				const fullCentroid = getCentroid(currentFrame.positions, new Set());

				// Compute the bounding sphere radius
				let maxDistance = 0;
				currentFrame.positions.forEach((position) => {
					maxDistance = Math.max(
						maxDistance,
						position.distanceTo(fullCentroid),
					);
				});

				// Assume a default FOV of 75 degrees for perspective camera
				const fov = (75 * Math.PI) / 180; // Convert FOV to radians
				let distance = maxDistance / Math.tan(fov / 2);
				// if distance is NaN or 0, use a default distance
				if (Number.isNaN(distance) || distance === 0) {
					distance = 10;
				}

				const resetCamera = {
					camera: new THREE.Vector3(distance, distance, distance),
					target: fullCentroid,
				};

				console.log("Resetting camera to original position:", resetCamera);
				setCameraAndControls(resetCamera);
			}
		};

		document.addEventListener("keydown", handleKeyDown);

		return () => {
			document.removeEventListener("keydown", handleKeyDown);
		};
	}, [
		setPlaying,
		playing,
		step,
		setStep,
		length,
		selectedFrames,
		bookmarks,
		setBookmarks,
		currentFrame,
		selectedIds,
		points,
		setPoints,
		selectedPoint,
		setSelectedPoint,
		showParticleInfo,
		setShowParticleInfo,
		isDrawing,
		setIsDrawing,
		cameraAndControls,
		setCameraAndControls,
	]);
};
