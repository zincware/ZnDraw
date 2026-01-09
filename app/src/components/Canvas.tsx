import * as THREE from "three";
import { useEffect, useRef, useState } from "react";
import { Canvas, useThree } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import type { OrbitControls as OrbitControlsImpl } from "three-stdlib";
import { useAppStore, getActiveCurves, selectPreferredCurve } from "../store";
import { useSettings } from "../hooks/useSettings";
import { useTheme } from "@mui/material/styles";
import {
	Snackbar,
	Alert,
	Box as MuiBox,
	CircularProgress,
	Typography,
	Button,
} from "@mui/material";
import {
	useCameraControls,
	type ControlsState,
} from "../hooks/useCameraControls";
import { useCameraSync, registerSession } from "../hooks/useCameraSync";
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
 * Handles three types of sync:
 * 1. Session camera sync (useCameraSync) - for Python client camera access
 * 2. Geometry camera sync (useGeometryCameraSync) - for updating Camera geometry from OrbitControls
 * 3. Geometry-to-controls sync - syncing OrbitControls.target when attached to geometry camera
 */
function CameraSyncIntegration({
	controlsRef,
	controlsState,
}: {
	controlsRef: React.RefObject<OrbitControlsImpl | null>;
	controlsState: ControlsState;
}) {
	const { camera } = useThree();
	const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);
	const geometries = useAppStore((state) => state.geometries);
	const curveRefs = useAppStore((state) => state.curveRefs);

	// Session camera sync (for Python client access)
	const { syncCamera } = useCameraSync(camera, controlsRef);

	// Geometry camera sync (for updating Camera geometry from OrbitControls)
	const { syncToGeometry, isUpdatingGeometry } = useGeometryCameraSync({
		camera,
		controlsRef,
		controlsState,
	});

	// Sync OrbitControls.target when attached to a geometry camera
	// This ensures OrbitControls knows the correct target when controls are enabled
	useEffect(() => {
		const controls = controlsRef.current;
		if (!controls) return;

		// Skip if we're currently pushing changes to geometry (prevent loops)
		if (isUpdatingGeometry.current) return;

		if (!attachedCameraKey) {
			// Not attached to any camera - OrbitControls.target stays as-is
			// This preserves the target from when we were attached
			return;
		}

		const cameraGeometry = geometries[attachedCameraKey];
		if (!cameraGeometry || cameraGeometry.type !== "Camera") return;

		const targetData = cameraGeometry.data?.target;
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
			} else {
				// Fallback: try to get single point from curve geometry
				const curveGeometry = geometries[curveKey];
				if (
					curveGeometry?.type === "Curve" &&
					curveGeometry.data?.position?.[0]
				) {
					targetPosition = curveGeometry.data.position[0];
				}
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
	}, [
		attachedCameraKey,
		geometries,
		curveRefs,
		controlsRef,
		isUpdatingGeometry,
	]);

	// Attach onChange to controls imperatively
	useEffect(() => {
		const controls = controlsRef.current;
		if (!controls) return;

		const handleChange = () => {
			// Sync to session camera (for Python access)
			syncCamera();
			// Sync to Camera geometry (if attached and editable)
			syncToGeometry();
		};

		controls.addEventListener("change", handleChange);
		return () => {
			controls.removeEventListener("change", handleChange);
		};
	}, [controlsRef, syncCamera, syncToGeometry]);

	return null;
}

// The main scene component
function MyScene() {
	// Use individual selectors to prevent unnecessary re-renders
	const roomId = useAppStore((state) => state.roomId);
	const sessionId = useAppStore((state) => state.sessionId);
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
	const mode = useAppStore((state) => state.mode);
	const theme = useTheme();

	// Track frame load time when not playing
	useFrameLoadTime();

	// Derive session camera key
	const sessionCameraKey = sessionId ? `cam:session:${sessionId}` : null;

	// Get the effective camera - attached camera or session camera
	const effectiveCameraKey = attachedCameraKey || sessionCameraKey;

	// Get camera control states based on effective camera
	const cameraControls = useCameraControls(effectiveCameraKey, geometries);

	// Ref for OrbitControls (used by camera sync)
	const orbitControlsRef = useRef<OrbitControlsImpl>(null);

	// Fetch settings (for lighting, pathtracing - NOT for camera anymore)
	const { data: settingsResponse, isLoading: settingsLoading } = useSettings(
		roomId || "",
	);

	// Connection state - needed to re-register session on socket reconnect
	const connected = useAppStore((state) => state.connected);

	// Register session on mount AND on socket reconnect
	// MUST be before loading check to avoid circular dependency
	// Session camera is created when session:register is called on backend
	// Re-registers when socket reconnects to ensure the new socket joins session room
	useEffect(() => {
		if (sessionId && connected) {
			const urlParams = new URLSearchParams(window.location.search);
			const alias = urlParams.get("alias");
			registerSession(sessionId, alias);
		}
	}, [sessionId, connected]);

	// Initialize attachedCameraKey to session camera when it becomes available
	// This ensures one camera is always active (radio button behavior)
	useEffect(() => {
		// Only set if not already set and session camera exists
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

	// Timeout state for session camera initialization
	const [sessionCameraTimeout, setSessionCameraTimeout] = useState(false);

	// Timeout effect: show error if session camera doesn't appear within 10 seconds
	useEffect(() => {
		// Reset timeout state when sessionId changes (e.g., reconnect)
		setSessionCameraTimeout(false);

		// Don't start timer if no sessionId or camera already exists
		if (!sessionId || sessionCameraData) return;

		const timer = setTimeout(() => {
			setSessionCameraTimeout(true);
		}, 10000); // 10 second timeout

		return () => clearTimeout(timer);
	}, [sessionId, sessionCameraData]);

	// Return early with loading state while settings or session camera are loading
	if (settingsLoading || !settingsResponse || !sessionCameraData) {
		// Show error state on timeout
		if (sessionCameraTimeout) {
			return (
				<MuiBox
					sx={{
						width: "100%",
						height: "calc(100vh - 64px)",
						display: "flex",
						flexDirection: "column",
						alignItems: "center",
						justifyContent: "center",
						bgcolor: theme.palette.background.default,
						gap: 2,
					}}
				>
					<Typography variant="h6" color="error">
						Failed to initialize camera session
					</Typography>
					<Typography variant="body2" color="text.secondary">
						The connection to the server may have been interrupted.
					</Typography>
					<Button
						variant="contained"
						onClick={() => window.location.reload()}
						sx={{ mt: 2 }}
					>
						Refresh Page
					</Button>
				</MuiBox>
			);
		}
		return (
			<MuiBox
				sx={{
					width: "100%",
					height: "calc(100vh - 64px)",
					display: "flex",
					alignItems: "center",
					justifyContent: "center",
					bgcolor: theme.palette.background.default,
				}}
			>
				<CircularProgress />
			</MuiBox>
		);
	}

	// Backend always returns defaults, so these are guaranteed to exist
	const studioLightingSettings = settingsResponse.data.studio_lighting;
	const pathtracingSettings = settingsResponse.data.pathtracing;
	const pathtracingEnabled = pathtracingSettings.enabled === true;

	// Get camera settings from session camera geometry
	const cameraPosition = sessionCameraData.position as [number, number, number];
	const cameraFov = sessionCameraData.fov ?? 75;
	const cameraType = sessionCameraData.camera_type ?? "PerspectiveCamera";
	const preserveDrawingBuffer =
		sessionCameraData.preserve_drawing_buffer ?? false;
	const showCrosshair = sessionCameraData.show_crosshair ?? false;

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
