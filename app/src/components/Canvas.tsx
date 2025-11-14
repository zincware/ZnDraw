import * as THREE from "three";
import { useEffect } from "react";
import { Canvas } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import { useAppStore, getActiveCurves, selectPreferredCurve } from "../store";
import { useExtensionData } from "../hooks/useSchemas";
import { useTheme } from "@mui/material/styles";
import { Snackbar, Alert, Box as MuiBox, CircularProgress } from "@mui/material";
import { useCameraControls } from "../hooks/useCameraControls";

// Import our new, self-contained components
import { Cell } from "./three/Cell";
import { Floor } from "./three/Floor";
import Sphere from "./three/Particles";
import Arrow from "./three/Arrow";
import Bonds from "./three/Bonds";
import Box from "./three/Box";
import Plane from "./three/Plane";
import Curve from "./three/Curve";
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

// The main scene component
function MyScene() {
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const userName = useAppStore((state) => state.userName);
  const geometries = useAppStore((state) => state.geometries);
  const activeCurveForDrawing = useAppStore((state) => state.activeCurveForDrawing);
  const setActiveCurveForDrawing = useAppStore((state) => state.setActiveCurveForDrawing);
  const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);
  const snackbar = useAppStore((state) => state.snackbar);
  const hideSnackbar = useAppStore((state) => state.hideSnackbar);
  const mode = useAppStore((state) => state.mode);
  const theme = useTheme();

  // Track frame load time when not playing
  useFrameLoadTime();

  // Get camera control states based on attached camera
  const cameraControls = useCameraControls(attachedCameraKey, geometries);

  const { data: studioLightingSettings } = useExtensionData(
    roomId || "",
    "settings",
    "studio_lighting",
  );

  const { data: cameraSettings } = useExtensionData(
    roomId || "",
    "settings",
    "camera",
  );

  const { data: pathtracingSettings } = useExtensionData(
    roomId || "",
    "settings",
    "pathtracing",
  );

  const pathtracingEnabled = pathtracingSettings?.enabled === true;

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

  // Return early with loading state if required settings are not yet loaded
  if (!studioLightingSettings || !cameraSettings) {
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

  const backgroundColor = studioLightingSettings.background_color === "default"
    ? theme.palette.background.default
    : studioLightingSettings.background_color;


  return (
    <div style={{ width: "100%", height: "calc(100vh - 64px)" }}>
      <Canvas
        // Add a key that changes when the camera type changes.
        // This will force React to re-create the Canvas and its camera.
        key={cameraSettings.camera}
        shadows
        // We REMOVE the dynamic near/far properties from here.
        // The initial position and fov are still useful.
        camera={{ position: [-10, 10, 30], fov: 50 }}
        gl={{
          antialias: true,
          toneMapping: THREE.ACESFilmicToneMapping,
          preserveDrawingBuffer: cameraSettings.preserve_drawing_buffer === true
        }}
        style={{ background: backgroundColor }}
        // The orthographic prop still sets the initial camera type.
        orthographic={cameraSettings.camera === "OrthographicCamera"}
      >
        {/* Place the CameraManager here, inside the Canvas */}
        <CameraManager settings={cameraSettings} />

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
                    <Curve
                      geometryKey={name}
                      data={config.data}
                    />
                  </GeometryErrorBoundary>
                );
              } else if (config.type === "Camera") {
                return (
                  <GeometryErrorBoundary key={name} geometryKey={name}>
                    <Camera
                      geometryKey={name}
                      data={config.data}
                    />
                  </GeometryErrorBoundary>
                );
              } else if (config.type === "Cell") {
                return (
                  <GeometryErrorBoundary key={name} geometryKey={name}>
                    <Cell
                      data={config.data}
                    />
                  </GeometryErrorBoundary>
                );
              } else if (config.type === "Floor") {
                return (
                  <GeometryErrorBoundary key={name} geometryKey={name}>
                    <Floor
                      data={config.data}
                    />
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
              } else {
                console.warn(`Unhandled geometry type: ${config.type}`);
                return null;
              }
            })}
          {cameraSettings.show_crosshair && (
            <Crosshair />
          )}
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
          enableDamping={false}
          makeDefault
          enabled={cameraControls.enabled}
          enablePan={cameraControls.enablePan}
          enableRotate={cameraControls.enableRotate}
          enableZoom={cameraControls.enableZoom}
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
