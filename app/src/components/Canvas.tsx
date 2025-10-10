import * as THREE from "three";
import { useEffect, useRef } from "react";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { OrbitControls, ContactShadows, Environment } from "@react-three/drei";
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
import CameraManager from "./CameraManager";
import SceneLighting from "./SceneLighting"

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
        // Add a key that changes when the camera type changes.
        // This will force React to re-create the Canvas and its camera.
        key={cameraSettings.camera}
        shadows
        // We REMOVE the dynamic near/far properties from here.
        // The initial position and fov are still useful.
        camera={{ position: [-10, 10, 30], fov: 50 }}
        gl={{ antialias: true, toneMapping: THREE.ACESFilmicToneMapping }}
        style={{ background: backgroundColor }}
        // The orthographic prop still sets the initial camera type.
        orthographic={cameraSettings.camera === "OrthographicCamera"}
      >
        {/* Place the CameraManager here, inside the Canvas */}
        <CameraManager settings={cameraSettings} />
        <SceneLighting 
          ambient_light={studioLightingSettings.ambient_light}
          key_light={studioLightingSettings.key_light}
          fill_light={studioLightingSettings.fill_light}
          rim_light={studioLightingSettings.rim_light}
          hemisphere_light={studioLightingSettings.hemisphere_light}
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

        {/* {studioLightingSettings.contact_shadow && (
          <mesh receiveShadow rotation={[-Math.PI / 2, 0, 0]} position={[0, -15, 0]}>
            <planeGeometry args={[1000, 1000]} />
            <shadowMaterial opacity={0.3} />
          </mesh>
        )} */}

        <OrbitControls enableDamping={false} makeDefault />
      </Canvas>
    </div>
  );
}

export default MyScene;
