import { TransformControls, Dodecahedron } from "@react-three/drei";
import { useEffect, useRef } from "react";
import * as THREE from "three";

export default function ControlsBuilder({
  points,
  setPoints,
  selectedPoint,
  setSelectedPoint,
}: {
  points: THREE.Vector3[];
  setPoints: any;
  selectedPoint: THREE.Vector3 | null;
  setSelectedPoint: any;
}) {
  const mesh = useRef<THREE.Object3D>(null); // TODO: check type
  const controls = useRef<THREE.TransformControls>(null);

  useEffect(() => {
    // event on pressing backspace
    const handleKeyDown = (event: KeyboardEvent) => {
      if (event.key === "Backspace") {
        if (selectedPoint !== null) {
          const newPoints = points.filter(
            (point) => point.distanceTo(selectedPoint) > 0.1,
          );
          setSelectedPoint(null);
          setPoints(newPoints);
        } else if (points.length > 0) {
          // pop last point from points
          setSelectedPoint(null);
          setPoints(points.slice(0, points.length - 1));
        }
      }
    };

    document.addEventListener("keydown", handleKeyDown);
    return () => {
      document.removeEventListener("keydown", handleKeyDown);
    };
  }, [points, selectedPoint]);

  useEffect(() => {
    if (selectedPoint == null) {
      return;
    }

    if (controls.current && mesh.current) {
      controls.current.attach(mesh.current);

      const handleChange = () => {
        if (mesh.current) {
          const index = points.findIndex(
            (point) => point.distanceTo(selectedPoint) < 0.1,
          );
          if (index == -1) {
            // TODO: check what would be best here?
            return;
          }

          const newPosition = mesh.current.position.clone();
          // update the position of points[index] to newPosition
          const newPoints = [...points];
          newPoints[index] = newPosition;
          setPoints(newPoints);
        }
      };

      controls.current.addEventListener("objectChange", handleChange);

      // Clean up event listeners on unmount
      return () => {
        if (controls.current) {
          controls.current.removeEventListener("objectChange", handleChange);
        }
      };
    }
  }, [selectedPoint]);

  return (
    <>
      {selectedPoint !== null && (
        <>
          <TransformControls ref={controls} />
          <Dodecahedron
            args={[0.01, 0]}
            position={selectedPoint}
            material-color="black"
            ref={mesh}
          />
        </>
      )}
    </>
  );
}
