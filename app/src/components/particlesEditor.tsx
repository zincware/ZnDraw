import { useEffect, useRef, useCallback, useMemo, useState } from "react";
import { TransformControls } from "@react-three/drei";
import { Vector3 } from "three";

export const getCentroid = (positions: Vector3[], selection: Set<number>) => {
  const centroid = new Vector3();
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
const useCentroid = ({ frame, selectedIds }: any) => {
  return useMemo(() => {
    return selectedIds.size > 0
      ? getCentroid(frame.positions, selectedIds)
      : new Vector3();
  }, [frame.positions, selectedIds]);
};

export const ParticleControls = ({
  frame,
  selectedIds,
  setFrame,
  highlight,
}) => {
  const controls = useRef(null);
  const controlsPostRef = useRef(new Vector3());

  // State for the edit mode: "None", "translate", or "rotate"
  const [mode, setMode] = useState("None");

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
          selectedIds.has(i) ? pos.clone().sub(deltaPosition) : pos
        ),
      }));
    },
    [setFrame, selectedIds]
  );

  // Handle control changes, applying only necessary updates to position and delta
  const handleControlsChange = useCallback(() => {
    if (controls.current?.object?.position && selectedIds.size > 0) {
      const deltaPosition = controlsPostRef.current
        .clone()
        .sub(controls.current.object.position);
      applyDeltaToPositions(deltaPosition);
      controlsPostRef.current.copy(controls.current.object.position);
    }
  }, [applyDeltaToPositions, selectedIds, mode]);

  // Toggle mode between "None", "translate", and "rotate" on "E" key press
  useEffect(() => {
    const toggleMode = (event) => {
      if (event.key.toLowerCase() === "e") {
        setMode((prevMode) => {
          switch (prevMode) {
            case "None":
              return "translate";
            case "translate":
              return "rotate";
            case "rotate":
              return "None";
            default:
              return "None";
          }
        });
      }
    };

    window.addEventListener("keydown", toggleMode);
    return () => {
      window.removeEventListener("keydown", toggleMode);
    };
  }, []);

  // Apply mode to TransformControls whenever it changes
  useEffect(() => {
    if (controls.current) {
      controls.current.mode = mode === "None" ? "" : mode;
    }
  }, [mode]);

  return (
    <>
      {highlight === "selection" && selectedIds.size > 0 && mode !== "None" && (
        <TransformControls ref={controls} onChange={handleControlsChange} />
      )}
    </>
  );
};