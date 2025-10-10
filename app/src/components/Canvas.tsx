import * as THREE from "three";
import { useEffect, useRef } from "react";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { OrbitControls, ContactShadows } from "@react-three/drei";
import { useAppStore } from "../store";
import { useExtensionData } from "../hooks/useSchemas";
import { useColorScheme, useTheme } from "@mui/material/styles";

// Import our new, self-contained components
import { Cell } from "./three/Cell";
import Sphere from "./three/Particles";
import Arrow from "./three/Arrow";
import Bonds from "./three/SingleBonds";
import Curve from "./three/Curve";
import VirtualCanvas from "./three/VirtualCanvas";
import Crosshair from "./three/crosshair";

// A small helper component to keep a light attached to the camera
function CameraAttachedLight({ intensity = 1.0 }) {
  const { camera } = useThree();
  const lightRef = useRef<THREE.DirectionalLight>(null!);

  useFrame(() => {
    if (lightRef.current) {
      lightRef.current.position.copy(camera.position);
    }
  });

  return <directionalLight ref={lightRef} intensity={intensity} />;
}

// The main scene component
function MyScene() {
  const { roomId, userId, geometries } = useAppStore();
  const { mode } = useColorScheme();

  const { data: studioLightingSettings } = useExtensionData(
    roomId || "",
    userId || "",
    "settings",
    "studio_lighting",
  );

  const { data: cameraSettings } = useExtensionData(
    roomId || "",
    userId || "",
    "settings",
    "camera",
  );

  useEffect(() => {
    console.log("camera-settings", cameraSettings);
  }, [cameraSettings])

  const backgroundColor = studioLightingSettings.background_color === "default" ? (mode === "light" ? "#FFFFFF" : "#212121ff") : studioLightingSettings.background_color;


  return (
    <div style={{ width: "100%", height: "calc(100vh - 64px)" }}>
      <Canvas
        shadows
        camera={{ position: [-10, 10, 30], fov: 50, near: cameraSettings.near_plane, far: cameraSettings.far_plane}}
        gl={{ antialias: true, toneMapping: THREE.ACESFilmicToneMapping }}
        style={{ background: backgroundColor }}
        orthographic={cameraSettings.camera === "OrthographicCamera"}
      >
        <CameraAttachedLight intensity={studioLightingSettings.key_light_intensity} />
        <ambientLight intensity={studioLightingSettings.fill_light_intensity} />
        <directionalLight
          position={[10, 20, -20]}
          intensity={studioLightingSettings.rim_light_intensity}
        />

        {/* Render our clean, refactored components */}
        {/* <Sphere /> */}
        {Object.entries(geometries)
          .filter(([_, config]) => config.data?.active !== false)
          .map(([name, config]) => {
            if (config.type === "Sphere") {
              return (
                <Sphere
                  key={name}
                  geometryKey={name}
                  data={config.data}
                />
              );
            } else if (config.type === "Bond") {
              return (
                <Bonds
                  key={name}
                  geometryKey={name}
                  data={config.data}
                />
              );
            } else if (config.type === "Arrow") {
              return (
                <Arrow
                  key={name}
                  geometryKey={name}
                  data={config.data}
                />
              );
            } else if (config.type === "Curve") {
              return (
                <Curve
                  key={name}
                  geometryKey={name}
                  data={config.data}
                />
              );
            } else if (config.type === "Cell") {
              return (
                <Cell
                  key={name}
                  data={config.data}
                />
              );
            } else {
              console.warn(`Unhandled geometry type: ${config.type}`);
              return null;
            }

            // You could add other types here
            // if (config.type === 'Box') {
            //   return <Box key={name} /* ...box props... */ />;
            // }

            return null;
          })}
        {cameraSettings.show_crosshair && (
          <Crosshair />
        )}
        <VirtualCanvas />

        {studioLightingSettings.contact_shadow && (
          <ContactShadows
            position={[0, -15, 0]}
            opacity={0.5}
            scale={80}
            blur={2}
            far={30}
            color="#000000"
          />
        )}

        <OrbitControls enableDamping={false} makeDefault />
      </Canvas>
    </div>
  );
}

export default MyScene;
