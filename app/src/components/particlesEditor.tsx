import { useEffect, useRef, useCallback, useMemo } from 'react';
import { TransformControls } from "@react-three/drei";
import { Vector3 } from 'three'; // assuming you're using three.js

export const getCentroid = (
    positions: Vector3[],
    selection: Set<number>,
  ) => {
    if (selection.size > 0) {
      const centroid = new Vector3();
      selection.forEach((i) => {
        centroid.add(positions[i]);
      });
      centroid.divideScalar(selection.size);
      return centroid;
    } else {
      const centroid = new Vector3();
      positions.forEach((position) => {
        centroid.add(position);
      });
      centroid.divideScalar(positions.length);
      return centroid;
    }
  };

// Custom hook for handling centroid calculations
const useCentroid = ({frame, selectedIds}: any) => {
  const centroid = useMemo(() => {
    return selectedIds.size > 0 ? getCentroid(frame.positions, selectedIds) : new Vector3();
  }, [frame.positions, selectedIds]);
  return centroid;
};

export const ParticleControls = ({ frame, selectedIds, setFrame, highlight }: any) => {
  const controls = useRef(null);
  const controlsPostRef = useRef(new Vector3());

  // Efficiently calculate centroid and attach control to it when `selectedIds` changes
  const centroid = useCentroid({frame, selectedIds});

  useEffect(() => {
    if (controls.current && selectedIds.size > 0) {
      controls.current.object.position.copy(centroid);
      controlsPostRef.current.copy(centroid);
    }
  }, [centroid]);

  // Helper to update frame positions based on delta
  const applyDeltaToPositions = useCallback((deltaPosition) => {
    setFrame((prevFrame) => ({
      ...prevFrame,
      positions: prevFrame.positions.map((pos, i) =>
        selectedIds.has(i) ? pos.clone().sub(deltaPosition) : pos
      ),
    }));
  }, [setFrame, selectedIds]);

  // Handle control changes, applying only necessary updates to position and delta
  const handleControlsChange = useCallback(() => {
    if (controls.current?.object?.position && selectedIds.size > 0) {
      const deltaPosition = controlsPostRef.current.clone().sub(controls.current.object.position);
      applyDeltaToPositions(deltaPosition);
      controlsPostRef.current.copy(controls.current.object.position);
    }
  }, [applyDeltaToPositions, selectedIds]);

  return (
    <>
      {highlight === 'selection' && selectedIds.size > 0 && (
        <TransformControls ref={controls} onChange={handleControlsChange} />
      )}
    </>
  );
};