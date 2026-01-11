import * as THREE from "three";
import { useEffect, useRef, type RefObject } from "react";
import { Canvas, useThree } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import type { OrbitControls as OrbitControlsImpl } from "three-stdlib";
import { useAppStore, getActiveCurves, selectPreferredCurve } from "../store";
import { useSettings } from "../hooks/useSettings";
import { useTheme } from "@mui/material/styles";
import { Snackbar, Alert } from "@mui/material";
import { CanvasLoadingState } from "./CanvasLoadingState";
import { CanvasErrorState } from "./CanvasErrorState";
import {
	useCameraControls,
	type ControlsState,
} from "../hooks/useCameraControls";
import { useGeometryCameraSync } from "../hooks/useGeometryCameraSync";

// Import our new, self-contained components
import { Cell } from "./three/Cell";
import { Floor } from "./three/Floor";
import Sphere from "./three/Particles";
import Arrow from "./three/Arrow";
import Bonds from "./three/Bonds";
import Box from "./three/Box";
import Plane from "./three/Plane";
import Curve from "./three/Curve";
import Shape from "./three/Shape";
import Camera from "./three/Camera";
import VirtualCanvas from "./three/VirtualCanvas";
import Crosshair from "./three/crosshair";
import CameraManager from "./CameraManager";
import SceneLighting from "./SceneLighting";
import { KeyboardShortcutsHandler } from "./three/KeyboardShortcutsHandler";
import StaticInfoBox from "./three/StaticInfoBox";
import HoverInfoBox from "./three/HoverInfoBox";
import DrawingIndicator from "./three/DrawingIndicator";
import EditingIndicator from "./three/EditingIndicator";
import MultiGeometryTransformControls from "./three/MultiGeometryTransformControls";
import { PathTracingRenderer } from "./PathTracingRenderer";
import { GeometryErrorBoundary } from "./three/GeometryErrorBoundary";
import { useFrameLoadTime } from "../hooks/useFrameLoadTime";

/**
 * Component for integrating camera sync inside the Canvas.
 * Must be a child of Canvas to use useThree.
 *
 * Handles two types of sync:
 * 1. Geometry camera sync (useGeometryCameraSync) - for updating Camera geometry from OrbitControls
 * 2. Geometry-to-controls sync - syncing OrbitControls.target when attached to geometry camera
 *
 * Note: Session cameras are now regular geometries. Python clients access them via
 * the geometry system, which broadcasts INVALIDATE_GEOMETRY to sync changes.
 */
function CameraSyncIntegration({
	controlsRef,
	controlsState,
}: {
	controlsRef: RefObject<OrbitControlsImpl | null>;
	controlsState: ControlsState;
}) {
	const { camera } = useThree();
	const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);
	const geometries = useAppStore((state) => state.geometries);
	const curveRefs = useAppStore((state) => state.curveRefs);

	// Geometry camera sync (for updating Camera geometry from OrbitControls)
	const { syncToGeometry, isEchoBack } = useGeometryCameraSync({
		camera,
		controlsRef,
		controlsState,
	});

	// Sync OrbitControls.target when attached to a geometry camera
	// This ensures OrbitControls knows the correct target when controls are enabled
	useEffect(() => {
		const controls = controlsRef.current;
		if (!controls) return;

		if (!attachedCameraKey) {
			return;
		}

		const cameraGeometry = geometries[attachedCameraKey];
		if (!cameraGeometry || cameraGeometry.type !== "Camera") return;

		const geomData = cameraGeometry.data;
		if (!geomData) return;

		// Skip if this is an echo-back of our own update (value-based detection)
		if (isEchoBack(geomData)) return;

		const targetData = geomData.target;
		if (!targetData) return;

		// Resolve target (either direct coords or CurveAttachment)
		let targetPosition: [number, number, number] | null = null;

		if (Array.isArray(targetData) && targetData.length === 3) {
			// Direct coordinates
			targetPosition = targetData as [number, number, number];
		} else if (
			targetData &&
			typeof targetData === "object" &&
			targetData.type === "curve_attachment"
		) {
			// CurveAttachment - resolve via curve
			const curveKey = targetData.geometry_key;
			const progress = targetData.progress || 0;
			const curve = curveRefs[curveKey];

			if (curve) {
				const point = curve.getPointAt(progress);
				targetPosition = [point.x, point.y, point.z];
			}
		}

		if (targetPosition) {
			controls.target.set(
				targetPosition[0],
				targetPosition[1],
				targetPosition[2],
			);
			controls.update();
		}
	}, [attachedCameraKey, geometries, curveRefs, controlsRef, isEchoBack]);

	// Attach onChange to controls imperatively
	useEffect(() => {
		const controls = controlsRef.current;
		if (!controls) return;

		const handleChange = () => {
			syncToGeometry();
		};

		controls.addEventListener("change", handleChange);
		return () => {
			controls.removeEventListener("change", handleChange);
		};
	}, [controlsRef, syncToGeometry]);

	return null;
}

// The main scene component
function MyScene() {
	const roomId = useAppStore((state) => state.roomId);
	const sessionId = useAppStore((state) => state.sessionId);
	const isConnected = useAppStore((state) => state.isConnected);
	const initializationError = useAppStore((state) => state.initializationError);
	const geometries = useAppStore((state) => state.geometries);
	const activeCurveForDrawing = useAppStore(
		(state) => state.activeCurveForDrawing,
	);
	const setActiveCurveForDrawing = useAppStore(
		(state) => state.setActiveCurveForDrawing,
	);
	const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);
	const attachToCamera = useAppStore((state) => state.attachToCamera);
	const snackbar = useAppStore((state) => state.snackbar);
	const hideSnackbar = useAppStore((state) => state.hideSnackbar);
	const theme = useTheme();

	useFrameLoadTime();

	const sessionCameraKey = sessionId ? `cam:session:${sessionId}` : null;

	const cameraControls = useCameraControls(attachedCameraKey, geometries);

	const orbitControlsRef = useRef<OrbitControlsImpl | null>(null);

	const { data: settingsResponse } = useSettings(roomId || "");
	useEffect(() => {
		if (
			!attachedCameraKey &&
			sessionCameraKey &&
			geometries[sessionCameraKey]
		) {
			attachToCamera(sessionCameraKey);
		}
	}, [attachedCameraKey, sessionCameraKey, geometries, attachToCamera]);

	// Auto-select default curve on startup
	useEffect(() => {
		// Only run if we don't have an active curve selected yet
		if (activeCurveForDrawing) return;

		// Get all active curves using helper function
		const activeCurves = getActiveCurves(geometries);

		// If no curves available, do nothing
		if (activeCurves.length === 0) return;

		// Select preferred curve using helper function
		const defaultCurve = selectPreferredCurve(activeCurves);

		setActiveCurveForDrawing(defaultCurve);
	}, [geometries, activeCurveForDrawing, setActiveCurveForDrawing]);

	// Get session camera geometry data
	const sessionCamera = sessionCameraKey ? geometries[sessionCameraKey] : null;
	const sessionCameraData = sessionCamera?.data;

	// Show error state if initialization failed
	if (initializationError) {
		return <CanvasErrorState error={initializationError} />;
	}

	// Return early with loading state until fully connected and data is ready
	// Gate on: 1) isConnected (socket connected), 2) sessionId (room joined),
	// 3) settingsResponse (settings loaded), 4) sessionCameraData (camera geometry loaded)
	if (!isConnected || !sessionId || !settingsResponse || !sessionCameraData) {
		return <CanvasLoadingState />;
	}

	// Backend always returns defaults, so these are guaranteed to exist
	const studioLightingSettings = settingsResponse.data.studio_lighting;
	const pathtracingSettings = settingsResponse.data.pathtracing;
	const pathtracingEnabled = pathtracingSettings.enabled === true;

	const cameraPosition = sessionCameraData.position as [number, number, number];
	const cameraFov = sessionCameraData.fov;
	const cameraType = sessionCameraData.camera_type;
	const preserveDrawingBuffer = sessionCameraData.preserve_drawing_buffer;
	const showCrosshair = sessionCameraData.show_crosshair;

	const backgroundColor =
		studioLightingSettings.background_color === "default"
			? theme.palette.background.default
			: studioLightingSettings.background_color;

	return (
		<div style={{ width: "100%", height: "calc(100vh - 64px)" }}>
			<Canvas
				// Add a key that changes when the camera type changes.
				// This will force React to re-create the Canvas and its camera.
				key={cameraType}
				shadows
				// Use session camera position from geometry
				camera={{
					position: cameraPosition,
					fov: cameraFov,
				}}
				gl={{
					antialias: true,
					toneMapping: THREE.ACESFilmicToneMapping,
					preserveDrawingBuffer: preserveDrawingBuffer,
				}}
				style={{ background: backgroundColor }}
				// The orthographic prop sets the initial camera type.
				orthographic={cameraType === "OrthographicCamera"}
			>
				{/* Place the CameraManager here, inside the Canvas */}
				<CameraManager sessionCameraData={sessionCameraData} />

				{/* Wrap scene in PathTracingRenderer */}
				<PathTracingRenderer settings={pathtracingSettings}>
					{/* Disable studio lighting when pathtracing (environment provides light) */}
					{!pathtracingEnabled && (
						<SceneLighting
							ambient_light={studioLightingSettings.ambient_light}
							key_light={studioLightingSettings.key_light}
							fill_light={studioLightingSettings.fill_light}
							rim_light={studioLightingSettings.rim_light}
							hemisphere_light={studioLightingSettings.hemisphere_light}
						/>
					)}

					{/* Keyboard shortcuts for 3D interactions */}
					<KeyboardShortcutsHandler />

					{/* Render our clean, refactored components */}
					{/* <Sphere /> */}
					{Object.entries(geometries)
						.filter(([_, config]) => config.data?.active !== false)
						.map(([name, config]) => {
							if (config.type === "Sphere") {
								return (
									<GeometryErrorBoundary key={name} geometryKey={name}>
										<Sphere
											geometryKey={name}
											data={config.data}
											pathtracingEnabled={pathtracingEnabled}
										/>
									</GeometryErrorBoundary>
								);
							} else if (config.type === "Bond") {
								return (
									<GeometryErrorBoundary key={name} geometryKey={name}>
										<Bonds
											geometryKey={name}
											data={config.data}
											pathtracingEnabled={pathtracingEnabled}
										/>
									</GeometryErrorBoundary>
								);
							} else if (config.type === "Arrow") {
								return (
									<GeometryErrorBoundary key={name} geometryKey={name}>
										<Arrow
											geometryKey={name}
											data={config.data}
											pathtracingEnabled={pathtracingEnabled}
										/>
									</GeometryErrorBoundary>
								);
							} else if (config.type === "Curve") {
								return (
									<GeometryErrorBoundary key={name} geometryKey={name}>
										<Curve geometryKey={name} data={config.data} />
									</GeometryErrorBoundary>
								);
							} else if (config.type === "Camera") {
								return (
									<GeometryErrorBoundary key={name} geometryKey={name}>
										<Camera geometryKey={name} data={config.data} />
									</GeometryErrorBoundary>
								);
							} else if (config.type === "Cell") {
								return (
									<GeometryErrorBoundary key={name} geometryKey={name}>
										<Cell data={config.data} />
									</GeometryErrorBoundary>
								);
							} else if (config.type === "Floor") {
								return (
									<GeometryErrorBoundary key={name} geometryKey={name}>
										<Floor data={config.data} />
									</GeometryErrorBoundary>
								);
							} else if (config.type === "Box") {
								return (
									<GeometryErrorBoundary key={name} geometryKey={name}>
										<Box
											geometryKey={name}
											data={config.data}
											pathtracingEnabled={pathtracingEnabled}
										/>
									</GeometryErrorBoundary>
								);
							} else if (config.type === "Plane") {
								return (
									<GeometryErrorBoundary key={name} geometryKey={name}>
										<Plane
											geometryKey={name}
											data={config.data}
											pathtracingEnabled={pathtracingEnabled}
										/>
									</GeometryErrorBoundary>
								);
							} else if (config.type === "Shape") {
								return (
									<GeometryErrorBoundary key={name} geometryKey={name}>
										<Shape
											geometryKey={name}
											data={config.data}
											pathtracingEnabled={pathtracingEnabled}
										/>
									</GeometryErrorBoundary>
								);
							} else {
								console.warn(`Unhandled geometry type: ${config.type}`);
								return null;
							}
						})}
					{showCrosshair && <Crosshair />}
					<VirtualCanvas />

					{/* Multi-geometry transform controls for editing mode */}
					<MultiGeometryTransformControls />
				</PathTracingRenderer>

				{/* {studioLightingSettings.contact_shadow && (
          <mesh receiveShadow rotation={[-Math.PI / 2, 0, 0]} position={[0, -15, 0]}>
            <planeGeometry args={[1000, 1000]} />
            <shadowMaterial opacity={0.3} />
          </mesh>
        )} */}

				<OrbitControls
					ref={orbitControlsRef}
					enableDamping={false}
					makeDefault
					enabled={cameraControls.enabled}
					enablePan={cameraControls.enablePan}
					enableRotate={cameraControls.enableRotate}
					enableZoom={cameraControls.enableZoom}
				/>

				{/* Camera sync integration for Python-side camera access and geometry sync */}
				<CameraSyncIntegration
					controlsRef={orbitControlsRef}
					controlsState={cameraControls}
				/>
			</Canvas>
			{/* Info boxes and drawing/editing indicators rendered outside Canvas, in DOM */}
			<StaticInfoBox />
			<HoverInfoBox />
			<DrawingIndicator />
			<EditingIndicator />

			{/* Global snackbar for notifications */}
			{snackbar && (
				<Snackbar
					open={snackbar.open}
					autoHideDuration={6000}
					onClose={hideSnackbar}
					anchorOrigin={{ vertical: "top", horizontal: "center" }}
				>
					<Alert severity={snackbar.severity} onClose={hideSnackbar}>
						{snackbar.message}
					</Alert>
				</Snackbar>
			)}
		</div>
	);
}

export default MyScene;
