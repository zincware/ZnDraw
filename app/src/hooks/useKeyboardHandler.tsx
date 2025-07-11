import { useEffect } from "react";
import { useAppContext } from "../contexts/AppContext";
import { getCentroid } from "../components/particlesEditor";

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

	useEffect(() => {
		const handleKeyDown = (event: KeyboardEvent) => {

			const canvasContainer = document.querySelector('.canvas-container');
			const isCanvasFocused = canvasContainer && canvasContainer.contains(document.activeElement);
			const isBodyFocused = document.activeElement === document.body;
			
			if (!isCanvasFocused && !isBodyFocused) {
				return;
			}

			if (event.key === " ") {
				console.log("Space bar pressed, current playing state:", playing);
				event.preventDefault();
				event.stopPropagation();
				setPlaying((prevPlaying) => {
					console.log("Toggling playing from", prevPlaying, "to", !prevPlaying);
					// Ensure the state is actually updated and returned
					const newState = !prevPlaying;
					console.log("New playing state after toggle:", newState);
					return newState;
				});
				return;
			}

			if (event.key === "ArrowRight") {
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
				// Delete selected atoms
				if (selectedIds.size > 0) {
					// Logic for deleting atoms would go here
					console.log("Delete selected atoms:", selectedIds);
				}
				// Delete selected point
				if (selectedPoint) {
					const updatedPoints = points.filter(
						(point) => point !== selectedPoint,
					);
					setPoints(updatedPoints);
					setSelectedPoint(null);
				}
			} else if (event.key === "Escape") {
				// Clear selection
				setSelectedPoint(null);
			} else if (event.key === "c" || event.key === "C") {
				// Center camera on centroid of selected particles, or all particles if none selected
				const centroid = selectedIds.size > 0 
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
