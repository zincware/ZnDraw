import { Dodecahedron, TransformControls } from "@react-three/drei";
import { useEffect, useRef } from "react";
import type * as THREE from "three";

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
	const controls = useRef<typeof TransformControls>(null);

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
