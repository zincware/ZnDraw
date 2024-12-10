import React from "react";
import { TransformControls } from "@react-three/drei";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { Euler, Vector3 } from "three";
import { socket } from "../socket";

export const getCentroid = (positions: Vector3[], selection: Set<number>) => {
	const centroid = new Vector3();
	if (!positions || positions.length === 0) {
		return centroid;
	}
	if (selection.size > 0) {
		selection.forEach((i) => {
			centroid.add(positions[i]);
		});
		centroid.divideScalar(selection.size);
	} else {
		positions.forEach((position) => {
			centroid.add(position);
		});
		centroid.divideScalar(positions.length);
	}
	return centroid;
};

// Custom hook for handling centroid calculations
export const useCentroid = ({ frame, selectedIds }: any) => {
	return useMemo(() => {
		return getCentroid(frame.positions, selectedIds);
	}, [frame.positions, selectedIds]);
};

interface ParticleControlsProps {
	frame: any;
	selectedIds: Set<number>;
	setFrame: (frame: any) => void;
	roomLock: boolean;
	editMode: string;
	setEditMode: (mode: string) => void;
}

export const ParticleControls: React.FC<ParticleControlsProps> = ({
	frame,
	selectedIds,
	setFrame,
	roomLock,
	editMode,
	setEditMode,
}: ParticleControlsProps) => {
	const controls = useRef(null);
	const controlsPostRef = useRef(new Vector3());
	const controlsRotationRef = useRef(new Vector3());

	// State for the edit mode: "None", "translate", or "rotate"

	useEffect(() => {
		if (roomLock) {
			setEditMode("none");
		}
	}, [roomLock]);

	// Efficiently calculate centroid and attach control to it when `selectedIds` changes
	const centroid = useCentroid({ frame, selectedIds });

	useEffect(() => {
		if (controls.current && selectedIds.size > 0) {
			controls.current.object.position.copy(centroid);
			controlsPostRef.current.copy(centroid);
		}
	}, [centroid]);

	// Helper to update frame positions based on delta
	const applyDeltaToPositions = useCallback(
		(deltaPosition) => {
			setFrame((prevFrame) => ({
				...prevFrame,
				positions: prevFrame.positions.map((pos, i) =>
					selectedIds.has(i) ? pos.clone().sub(deltaPosition) : pos,
				),
			}));
		},
		[setFrame, selectedIds],
	);

	// Helper to update frame rotations based on delta rotation
	const applyDeltaToRotations = useCallback(
		(deltaRotation) => {
			setFrame((prevFrame) => ({
				...prevFrame,
				positions: prevFrame.positions.map((rot, i) => {
					if (selectedIds.has(i)) {
						// rotate the position around the centroid
						const position = rot.clone().sub(centroid);
						const euler = new Euler().setFromVector3(deltaRotation);
						position.applyEuler(euler);
						position.add(centroid);
						return position;
					}
					return rot;
				}),
			}));
		},
		[setFrame, selectedIds, centroid],
	);

	// Handle control changes, applying only necessary updates to position and delta
	const handleControlsChange = useCallback(() => {
		if (editMode === "translate") {
			if (controls.current?.object?.position && selectedIds.size > 0) {
				const deltaPosition = controlsPostRef.current
					.clone()
					.sub(controls.current.object.position);
				applyDeltaToPositions(deltaPosition);
				controlsPostRef.current.copy(controls.current.object.position);
			}
		} else if (editMode === "rotate") {
			if (controls.current?.object?.rotation && selectedIds.size > 0) {
				const deltaRotation = controlsRotationRef.current
					.clone()
					.sub(controls.current.object.rotation);
				applyDeltaToRotations(deltaRotation);
				controlsRotationRef.current.copy(controls.current.object.rotation);
			}
		}
	}, [applyDeltaToPositions, selectedIds, editMode, applyDeltaToRotations]);

	// Toggle mode between "None", "translate", and "rotate" on "E" key press
	useEffect(() => {
		const toggleMode = (event) => {
			if (document.activeElement !== document.body) {
				return;
			}
			if (roomLock) {
				return;
			}
			if (event.key.toLowerCase() === "e") {
				socket.emit("room:copy");
				setEditMode((prevMode) => {
					switch (prevMode) {
						case "none":
							return "translate";
						case "translate":
							return "rotate";
						case "rotate":
							return "none";
						default:
							return "none";
					}
				});
			}
		};

		window.addEventListener("keydown", toggleMode);
		return () => {
			window.removeEventListener("keydown", toggleMode);
		};
	}, [roomLock]);

	// Apply mode to TransformControls whenever it changes
	useEffect(() => {
		if (controls.current) {
			controls.current.mode = editMode === "none" ? "" : editMode;
		}
	}, [editMode]);

	return (
		<>
			{selectedIds.size > 0 && editMode !== "none" && (
				<TransformControls ref={controls} onChange={handleControlsChange} />
			)}
		</>
	);
};
