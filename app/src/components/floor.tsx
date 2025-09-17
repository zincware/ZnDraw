import { Plane } from "@react-three/drei";
import { Line } from "@react-three/drei";
import type React from "react";
import { useEffect, useState } from "react";
import * as THREE from "three";

function Grid({
	position = [0, 0, 0],
	gridSpacing = 10,
	sizeX = 500,
	sizeY = 500,
	color = "black",
}) {
	const lines = [];

	// Generate vertical lines
	for (let x = -sizeX / 2; x <= sizeX / 2; x += gridSpacing) {
		lines.push(
			<Line
				key={`v-${x}`}
				points={[
					new THREE.Vector3(x, 0, -sizeY / 2),
					new THREE.Vector3(x, 0, sizeY / 2),
				]}
				color={color}
				fog
			/>,
		);
	}

	// Generate horizontal lines
	for (let y = -sizeY / 2; y <= sizeY / 2; y += gridSpacing) {
		lines.push(
			<Line
				key={`h-${y}`}
				points={[
					new THREE.Vector3(-sizeX / 2, 0, y),
					new THREE.Vector3(sizeX / 2, 0, y),
				]}
				color={color}
				fog
			/>,
		);
	}

	return <group position={position}>{lines}</group>;
}

interface FloorProps {
	colorMode: string;
	roomConfig: {
		Visualization?: {
			simulation_box?: boolean;
			[key: string]: unknown;
		};
		[key: string]: unknown;
	};
}

export const Floor: React.FC<FloorProps> = ({ colorMode, roomConfig }) => {
	const [bsColor, setBsColor] = useState({
		"--bs-body-bg": "#fff",
		"--bs-secondary": "#fff",
	});

	useEffect(() => {
		// Re-read CSS custom properties when color mode changes
		console.log("Color mode changed to:", colorMode);
		const computedStyle = getComputedStyle(document.documentElement);
		setBsColor({
			"--bs-body-bg": computedStyle.getPropertyValue("--bs-body-bg"),
			"--bs-secondary": computedStyle.getPropertyValue("--bs-secondary"),
		});
	}, [colorMode]);

	return (
		<>
			{" "}
			<fog
				attach="fog"
				color={bsColor["--bs-body-bg"]}
				near={roomConfig.Camera.camera_far * 0.8}
				far={roomConfig.Camera.camera_far}
			/>
			<Plane
				args={[
					roomConfig.Camera.camera_far * 2,
					roomConfig.Camera.camera_far * 2,
					1,
					1,
				]}
				rotation={[-Math.PI / 2, 0, 0]}
				position={[0, -5, 0]}
				receiveShadow
			>
				<meshStandardMaterial color={bsColor["--bs-body-bg"]} />
			</Plane>
			<Grid
				position={[0, -4.95, 0]}
				gridSpacing={10}
				sizeX={roomConfig.Camera.camera_far * 2}
				sizeY={roomConfig.Camera.camera_far * 2}
				color={bsColor["--bs-secondary"]}
			/>
		</>
	);
};
