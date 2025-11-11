import { useEffect, useRef } from "react";
import { ContactShadows } from "@react-three/drei";
import { useThree } from "@react-three/fiber";
import * as THREE from "three";
import { useColorScheme, useTheme } from "@mui/material/styles";
import { useExtensionData } from "../../hooks/useSchemas";
import { useAppStore } from "../../store";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

interface FloorData {
  active: boolean;
  position: [number, number, number];
  color: string;
  height: number;
  grid_spacing: number;
  grid_color: string;
  grid_opacity: number;
  show_grid: boolean;
  show_shadows: boolean;
  shadow_opacity: number;
  shadow_blur: number;
  fog_enabled: boolean;
}

export const Floor = ({ data }: { data: FloorData }) => {
  const { scene } = useThree();
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const userName = useAppStore((state) => state.userName);
  const geometryDefaults = useAppStore((state) => state.geometryDefaults);
  const gridRef = useRef<THREE.GridHelper | null>(null);
  const { mode } = useColorScheme();
  const theme = useTheme();

  // Merge with defaults from Pydantic (single source of truth)
  const fullData = getGeometryWithDefaults<FloorData>(data, "Floor", geometryDefaults);

  // Get camera settings for camera_far
  const { data: cameraSettings } = useExtensionData(
    roomId || "",
    "settings",
    "camera"
  );

  // Get background color for fog
  const { data: studioSettings } = useExtensionData(
    roomId || "",
    "settings",
    "studio_lighting"
  );

  // Calculate floor size from camera far plane (use 2x for full coverage)
  const size = (cameraSettings?.far_plane || 300) * 2;

  // Cap shadow scale to maintain quality regardless of floor size
  // Shadows become too sparse if scale is too large relative to resolution
  const shadowScale = Math.min(size * 0.3, 150);

  // Map "default" colors to theme colors with better contrast
  const floorColor = fullData.color === "default"
    ? (mode === "light" ? "#e0e0e0" : "#303030")
    : fullData.color;

  const gridColor = fullData.grid_color === "default"
    ? (mode === "light" ? "#616161" : "#9e9e9e")
    : fullData.grid_color;

  // Sync fog with camera_far and background
  useEffect(() => {
    if (!fullData.fog_enabled) {
      scene.fog = null;
      return;
    }

    const bgColor = studioSettings?.background_color === "default"
      ? (mode === "light" ? "#FFFFFF" : "#212121")
      : studioSettings?.background_color || "#FFFFFF";
    const farPlane = cameraSettings?.far_plane || 300;

    // Fog extends from 60% to 100% of far plane
    const fogNear = farPlane * 0.6;
    const fogFar = farPlane;

    scene.fog = new THREE.Fog(bgColor, fogNear, fogFar);

    return () => {
      scene.fog = null;
    };
  }, [
    fullData.fog_enabled,
    cameraSettings?.far_plane,
    studioSettings?.background_color,
    mode,
    scene
  ]);

  // Update grid when settings change
  useEffect(() => {
    if (!gridRef.current || !fullData.show_grid) return;

    const grid = gridRef.current;
    grid.position.y = fullData.height;

    // Update grid material color/opacity
    if (grid.material instanceof THREE.Material) {
      grid.material.opacity = fullData.grid_opacity;
      grid.material.transparent = true;
      grid.material.color.set(gridColor);
    }
  }, [fullData.height, gridColor, fullData.grid_opacity, fullData.show_grid]);

  const divisions = Math.floor(size / fullData.grid_spacing);

  return (
    <group>
      {/* Floor Plane (receives shadows) */}
      <mesh
        rotation={[-Math.PI / 2, 0, 0]}
        position={[0, fullData.height - 0.01, 0]}
        receiveShadow
      >
        <planeGeometry args={[size, size]} />
        <meshStandardMaterial
          color={floorColor}
          roughness={0.8}
          metalness={0.2}
        />
      </mesh>

      {/* Grid Helper */}
      {fullData.show_grid && (
        <gridHelper
          ref={gridRef}
          args={[size, divisions, gridColor, gridColor]}
          position={[0, fullData.height, 0]}
        />
      )}

      {/* Contact Shadows (drei component) */}
      {fullData.show_shadows && (
        <ContactShadows
          position={[0, fullData.height + 0.01, 0]}
          opacity={fullData.shadow_opacity}
          scale={shadowScale}
          blur={fullData.shadow_blur}
          far={20}
          resolution={512}
        />
      )}
    </group>
  );
};
