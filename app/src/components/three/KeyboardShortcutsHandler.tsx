import { useEffect } from "react";
import { useThree } from "@react-three/fiber";
import { useQueryClient } from "@tanstack/react-query";
import { useAppStore } from "../../store";
import * as THREE from "three";
import { TypedArray } from "../../myapi/client";
import { OrbitControls as OrbitControlsImpl } from "three-stdlib";

// for now, we fix the key for positions
const positionsKey = `arrays.positions`;


const computeCentroid = (positions: TypedArray) => {
  const centroid = new THREE.Vector3(0, 0, 0);
  const count = positions.length / 3; // Assuming 3 components per vertex

  for (let i = 0; i < positions.length; i += 3) {
    centroid.x += positions[i];
    centroid.y += positions[i + 1];
    centroid.z += positions[i + 2];
  }

  centroid.divideScalar(count);
  return centroid;
};

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
    updateSelectionForGeometry,
    setIsDrawing,
    isDrawing,
    roomId,
    toggleInfoBoxes,
    geometries,
    activeCurveForDrawing,
    setActiveCurveForDrawing,
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
        let positionsData = queryClient.getQueryData<TypedArray | { [positionsKey]: TypedArray }>([
          "frame",
          roomId,
          currentFrame,
          positionsKey,
        ]);
        let positions: TypedArray | undefined;
        if (positionsData && typeof positionsData === "object" && positionsKey in positionsData) {
          positions = (positionsData as { [positionsKey]: TypedArray })[positionsKey];
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
        let positionsData = queryClient.getQueryData<TypedArray | { [positionsKey]: TypedArray }>([
          "frame",
          roomId,
          currentFrame,
          positionsKey,
        ]);
        let positions: TypedArray | undefined;
        if (positionsData && typeof positionsData === "object" && positionsKey in positionsData) {
          positions = (positionsData as { [positionsKey]: TypedArray })[positionsKey];
        } else if (positionsData) {
          positions = positionsData as TypedArray;
        }
        // filter positions by selection, if len(selection) > 0
        if (positions && selections["particles"] && selections["particles"].length > 0) {
          const filtered = new (positions.constructor as any)(selections["particles"].length * 3);
          selections["particles"].forEach((idx, i) => {
            if (positions) {
              filtered[i * 3] = positions[idx * 3];
              filtered[i * 3 + 1] = positions[idx * 3 + 1];
              filtered[i * 3 + 2] = positions[idx * 3 + 2];
            }
          });
          positions = filtered;
        } 
        if (positions && positions.length >= 3) {
          const centroid = computeCentroid(positions);
          if (controls) {
            const camera = controls.object;
            const target = controls.target;

            // Move camera to look at centroid from a distance
            const direction = new THREE.Vector3();
            camera.getWorldDirection(direction);
            const distance = camera.position.distanceTo(target);
            const newPosition = centroid.clone().add(direction.multiplyScalar(-distance));

            camera.position.copy(newPosition);
            controls.target.copy(centroid);
            controls.update();
          } else {
            console.warn("Camera controls not available");
          }
        }
        return;
      }

      // Handle x/X for toggle drawing with smart curve selection
      if (event.key === "x" || event.key === "X") {
        event.preventDefault();

        // Get all active curves
        const activeCurves = Object.entries(geometries)
          .filter(([_, g]) => g.type === "Curve" && g.data?.active !== false)
          .map(([key, _]) => key);

        if (activeCurves.length === 0) {
          // No curves available - do nothing or show notification
          console.warn("No active curves available for drawing");
          return;
        }

        if (activeCurves.length === 1) {
          // Only one curve - auto-select it
          const singleCurveKey = activeCurves[0];
          if (!activeCurveForDrawing) {
            setActiveCurveForDrawing(singleCurveKey);
          }
          setIsDrawing(!isDrawing);
          return;
        }

        // Multiple curves exist
        if (activeCurveForDrawing) {
          // Already have a selected curve - just toggle drawing
          setIsDrawing(!isDrawing);
        } else {
          // No curve selected - auto-select the first one (most recently added is last, but first is fine)
          setActiveCurveForDrawing(activeCurves[0]);
          setIsDrawing(true);
        }
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

    window.addEventListener("keydown", handleKeyDown);
    return () => window.removeEventListener("keydown", handleKeyDown);
  }, [
    controls,
    selections,
    updateSelectionForGeometry,
    isDrawing,
    setIsDrawing,
    queryClient,
    roomId,
    currentFrame,
    toggleInfoBoxes,
    geometries,
    activeCurveForDrawing,
    setActiveCurveForDrawing,
  ]);

  return null; // This component only handles keyboard events
};
