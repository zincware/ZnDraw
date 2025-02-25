import {
	OrbitControls,
	OrthographicCamera,
	PerspectiveCamera,
	TrackballControls,
} from "@react-three/drei";
import { Box } from "@react-three/drei";
import { debounce } from "lodash";
import {
	forwardRef,
	useCallback,
	useEffect,
	useMemo,
	useRef,
	useState,
} from "react";
import * as THREE from "three";
import { getCentroid, useCentroid } from "./particlesEditor";

const zeroVector = new THREE.Vector3(0, 0, 0);
const initialCameraVector = new THREE.Vector3(5, 5, 5);
const upVector = new THREE.Vector3(0, 1, 0);

const MoveCameraTarget = forwardRef(
	({ colorMode }: { colorMode: string }, ref: React.Ref<THREE.Group>) => {
		const shortDimension = 0.05;
		const longDimension = 0.5;

		return (
			<group ref={ref}>
				{/* X axis box */}
				<Box scale={[longDimension, shortDimension, shortDimension]}>
					<meshStandardMaterial
						color={colorMode === "light" ? "#454b66" : "#f5fdc6"}
					/>
				</Box>

				{/* Y axis box */}
				<Box scale={[shortDimension, longDimension, shortDimension]}>
					<meshStandardMaterial
						color={colorMode === "light" ? "#454b66" : "#f5fdc6"}
					/>
				</Box>

				{/* Z axis box */}
				<Box scale={[shortDimension, shortDimension, longDimension]}>
					<meshStandardMaterial
						color={colorMode === "light" ? "#454b66" : "#f5fdc6"}
					/>
				</Box>
			</group>
		);
	},
);

type CameraAndControls = {
	camera: THREE.Vector3;
	target: THREE.Vector3;
};

type CameraAndControlsProps = {
	cameraConfig: any;
	cameraAndControls: CameraAndControls;
	setCameraAndControls: React.Dispatch<React.SetStateAction<CameraAndControls>>;
	currentFrame: any;
	selectedIds: Set<number>;
	colorMode: string;
};

const CameraAndControls: React.FC<CameraAndControlsProps> = ({
	cameraConfig,
	cameraAndControls,
	setCameraAndControls,
	currentFrame,
	selectedIds,
	colorMode,
}) => {
	const cameraRef = useRef<any>(null);
	const controlsRef = useRef<any>(null);
	const centroid = useCentroid({
		frame: currentFrame,
		selectedIds: selectedIds,
	});

	// need this extra for the crosshair
	const crossHairRef = useRef<any>(null);

	const controlsOnChangeFn = useCallback((e: any) => {
		if (!crossHairRef.current) {
			return;
		}
		crossHairRef.current.position.copy(e.target.target);
	}, []);

	const controlsOnEndFn = useCallback(
		debounce(() => {
			if (!cameraRef.current || !controlsRef.current) {
				return;
			}
			setCameraAndControls({
				camera: cameraRef.current.position,
				target: controlsRef.current.target,
			});
		}, 100),
		[],
	);

	const rollCamera = useCallback((angle: number) => {
		if (!cameraRef.current) {
			return;
		}
		// Get the current direction the camera is looking
		const looksTo = new THREE.Vector3();
		cameraRef.current.getWorldDirection(looksTo);

		// Calculate the rotation axis (perpendicular to both yDir and looksTo)
		const rotationAxis = new THREE.Vector3();
		rotationAxis.crossVectors(upVector, looksTo).normalize();

		// Compute the quaternion for the roll rotation
		const quaternion = new THREE.Quaternion();
		quaternion.setFromAxisAngle(looksTo, angle);

		// Update the `up` vector of the camera
		const newUp = new THREE.Vector3();
		newUp.copy(cameraRef.current.up).applyQuaternion(quaternion).normalize();
		cameraRef.current.up.copy(newUp);

		// Update the camera matrix and controls
		cameraRef.current.updateProjectionMatrix();
		if (controlsRef.current) {
			controlsRef.current.update();
		}
	}, []);

	const getResetCamera = useCallback(() => {
		if (currentFrame.positions.length === 0) {
			return;
		}
		if (cameraRef.current === null) {
			return;
		}
		// now calculate the camera positions
		const fullCentroid = getCentroid(currentFrame.positions, new Set());

		// Compute the bounding sphere radius
		let maxDistance = 0;
		currentFrame.positions.forEach((x) => {
			maxDistance = Math.max(maxDistance, x.distanceTo(fullCentroid));
		});

		const fov = (cameraRef.current.fov * Math.PI) / 180; // Convert FOV to radians
		let distance = maxDistance / Math.tan(fov / 2);
		// if distance is NaN, return - happens for OrthographicCamera
		if (Number.isNaN(distance)) {
			distance = 0;
		}

		return {
			camera: new THREE.Vector3(distance, distance, distance),
			target: fullCentroid,
		};
	}, [currentFrame.positions]);

	// if the camera positions and target positions is default, adapt them to the scene
	useEffect(() => {
		if (
			cameraAndControls.camera.equals(initialCameraVector) &&
			cameraAndControls.target.equals(zeroVector)
		) {
			const resetCamera = getResetCamera();
			if (resetCamera) {
				setCameraAndControls(resetCamera);
			}
		}
	}, [currentFrame.positions]);

	// if the camera changes, run resetCamera
	useEffect(() => {
		const resetCamera = getResetCamera();
		if (resetCamera) {
			setCameraAndControls(resetCamera);
		}
	}, [cameraConfig.camera]);

	useEffect(() => {
		if (!cameraRef.current || !controlsRef.current) {
			return;
		}
		cameraRef.current.position.copy(cameraAndControls.camera);
		controlsRef.current.target.copy(cameraAndControls.target);
		controlsRef.current.update();
	}, [cameraAndControls]);

	// keyboard controls
	useEffect(() => {
		// page initialization

		const handleKeyDown = (event: KeyboardEvent) => {
			// if canvas is not focused, don't do anything
			if (document.activeElement !== document.body) {
				return;
			}
			if (event.key === "c") {
				setCameraAndControls((prev: any) => ({
					...prev,
					target: centroid,
				}));
			} else if (event.key === "o") {
				const resetCamera = getResetCamera();
				if (resetCamera) {
					setCameraAndControls(resetCamera);
				}
				// reset the camera roll
				if (cameraRef.current) {
					cameraRef.current.up.copy(upVector);
					cameraRef.current.updateProjectionMatrix();
					if (controlsRef.current) {
						controlsRef.current.update();
					}
				}
			} else if (event.key === "r") {
				const roll = Math.PI / 100;
				if (event.ctrlKey) {
					rollCamera(-roll);
				} else {
					rollCamera(roll);
				}
			}
		};

		// Add the event listener
		window.addEventListener("keydown", handleKeyDown);

		// Clean up the event listener on unmount
		return () => {
			window.removeEventListener("keydown", handleKeyDown);
		};
	}, [currentFrame, selectedIds]);

	return (
		<>
			{cameraConfig.camera === "OrthographicCamera" && (
				<OrthographicCamera
					makeDefault
					ref={cameraRef}
					zoom={10}
					near={cameraConfig.camera_near}
					far={cameraConfig.camera_far}
				>
					<pointLight decay={0} intensity={Math.PI / 2} />
				</OrthographicCamera>
			)}
			{cameraConfig.camera === "PerspectiveCamera" && (
				<PerspectiveCamera
					makeDefault
					ref={cameraRef}
					near={cameraConfig.camera_near}
					far={cameraConfig.camera_far}
				>
					<pointLight decay={0} intensity={Math.PI / 2} />
				</PerspectiveCamera>
			)}
			{cameraConfig.controls === "OrbitControls" && (
				<OrbitControls
					makeDefault
					onEnd={controlsOnEndFn}
					onChange={controlsOnChangeFn}
					enableDamping={false}
					ref={controlsRef}
				/>
			)}
			{cameraConfig.controls === "TrackballControls" && (
				<TrackballControls
					makeDefault
					onEnd={controlsOnEndFn}
					onChange={controlsOnChangeFn}
					ref={controlsRef}
				/>
			)}
			{cameraConfig.crosshair && controlsRef.current.target && (
				<MoveCameraTarget colorMode={colorMode} ref={crossHairRef} />
			)}
		</>
	);
};

export default CameraAndControls;
