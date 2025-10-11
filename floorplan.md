
  ---
  🏗️ Architecture Overview

  Backend (Python/Pydantic)          Frontend (React/Three.js)
  ├── src/zndraw/geometries/        ├── app/src/components/three/
  │   └── floor.py                   │   ├── Floor.tsx
  ├── Integration:                   ├── Integration:
  │   ├── __init__.py (registry)     │   └── Canvas.tsx (render)
  │   └── routes.py (default)        └── Fog + Shadows

  ---
  📝 Part 1: Backend - Pydantic Model

  File: src/zndraw/geometries/floor.py

  from pydantic import Field
  from .base import BaseGeometry

  class Floor(BaseGeometry):
      """A floor plane with grid, shadows, and fog support."""

      # Override base properties - Floor doesn't need dynamic data
      position: tuple[float, float, float] = Field(
          default=(0, 0, 0),
          description="Floor center position [x,y,z]"
      )

      color: str = Field(
          default="#808080",
          description="Floor base color (hex)"
      )

      # Floor-specific properties
      height: float = Field(
          default=0.0,
          description="Y-position of the floor plane"
      )

      size: float = Field(
          default=100.0,
          ge=10.0,
          le=1000.0,
          description="Floor size (square dimension)"
      )

      grid_spacing: float = Field(
          default=1.0,
          ge=0.5,
          le=50.0,
          description="Spacing between grid lines"
      )

      grid_color: str = Field(
          default="#404040",
          description="Grid line color (hex)"
      )

      grid_opacity: float = Field(
          default=0.5,
          ge=0.0,
          le=1.0,
          description="Grid line opacity"
      )

      show_grid: bool = Field(
          default=True,
          description="Toggle grid visibility"
      )

      show_shadows: bool = Field(
          default=True,
          description="Toggle contact shadows"
      )

      shadow_opacity: float = Field(
          default=0.5,
          ge=0.0,
          le=1.0,
          description="Shadow opacity"
      )

      shadow_blur: float = Field(
          default=2.0,
          ge=0.5,
          le=10.0,
          description="Shadow blur radius"
      )

      fog_enabled: bool = Field(
          default=True,
          description="Enable fog (synced with camera_far)"
      )

      fog_density: float = Field(
          default=0.5,
          ge=0.0,
          le=1.0,
          description="Fog density factor"
      )

  ---
  📝 Part 2: Frontend - Three.js Component

  File: app/src/components/three/Floor.tsx

  import { useEffect, useRef } from "react";
  import { ContactShadows } from "@react-three/drei";
  import { useThree } from "@react-three/fiber";
  import * as THREE from "three";
  import { useExtensionData } from "../../hooks/useSchemas";
  import { useAppStore } from "../../store";

  interface FloorData {
    active: boolean;
    position: [number, number, number];
    color: string;
    height: number;
    size: float;
    grid_spacing: number;
    grid_color: string;
    grid_opacity: number;
    show_grid: boolean;
    show_shadows: boolean;
    shadow_opacity: number;
    shadow_blur: number;
    fog_enabled: boolean;
    fog_density: number;
  }

  export const Floor = ({ data }: { data: FloorData }) => {
    const { scene } = useThree();
    const { roomId, userId } = useAppStore();
    const gridRef = useRef<THREE.GridHelper | null>(null);

    // Get camera settings for camera_far
    const { data: cameraSettings } = useExtensionData(
      roomId || "",
      userId || "",
      "settings",
      "camera"
    );

    // Get background color for fog
    const { data: studioSettings } = useExtensionData(
      roomId || "",
      userId || "",
      "settings",
      "studio_lighting"
    );

    // Sync fog with camera_far and background
    useEffect(() => {
      if (!data.fog_enabled) {
        scene.fog = null;
        return;
      }

      const bgColor = studioSettings?.background_color || "#FFFFFF";
      const farPlane = cameraSettings?.far_plane || 1000;

      // Fog distance controlled by fog_density (0-1)
      // Maps density: 0.0 = fog at far plane, 1.0 = fog at near
      const fogNear = farPlane * (1 - data.fog_density) * 0.5;
      const fogFar = farPlane;

      scene.fog = new THREE.Fog(bgColor, fogNear, fogFar);

      return () => {
        scene.fog = null;
      };
    }, [
      data.fog_enabled,
      data.fog_density,
      cameraSettings?.far_plane,
      studioSettings?.background_color,
      scene
    ]);

    // Update grid when settings change
    useEffect(() => {
      if (!gridRef.current || !data.show_grid) return;

      const grid = gridRef.current;
      grid.position.y = data.height;

      // Update grid material color/opacity
      if (grid.material instanceof THREE.Material) {
        grid.material.opacity = data.grid_opacity;
        grid.material.transparent = true;
        grid.material.color.set(data.grid_color);
      }
    }, [data.height, data.grid_color, data.grid_opacity, data.show_grid]);

    const divisions = Math.floor(data.size / data.grid_spacing);

    return (
      <group>
        {/* Floor Plane (receives shadows) */}
        <mesh 
          rotation={[-Math.PI / 2, 0, 0]} 
          position={[0, data.height - 0.01, 0]}
          receiveShadow
        >
          <planeGeometry args={[data.size, data.size]} />
          <meshStandardMaterial 
            color={data.color}
            roughness={0.8}
            metalness={0.2}
          />
        </mesh>

        {/* Grid Helper */}
        {data.show_grid && (
          <gridHelper
            ref={gridRef}
            args={[data.size, divisions, data.grid_color, data.grid_color]}
            position={[0, data.height, 0]}
          />
        )}

        {/* Contact Shadows (drei component) */}
        {data.show_shadows && (
          <ContactShadows
            position={[0, data.height + 0.01, 0]}
            opacity={data.shadow_opacity}
            scale={data.size}
            blur={data.shadow_blur}
            far={20}
            resolution={512}
          />
        )}
      </group>
    );
  };

  ---
  📝 Part 3: Integration Steps

  3.1. Register Floor Geometry

  File: src/zndraw/geometries/__init__.py

  from zndraw.geometries.floor import Floor

  geometries = {
      "Sphere": Sphere,
      "Arrow": Arrow,
      "Bond": Bond,
      "Curve": Curve,
      "Cell": Cell,
      "Floor": Floor,  # ← Add this
  }

  __all__ = [..., "Floor"]  # ← Add to exports

  3.2. Add Default Floor in /join

  File: src/zndraw/app/routes.py (line ~2304)

  from zndraw.geometries import Sphere, Bond, Curve, Cell, Floor  # ← Add Floor

  # After creating default geometries:
  r.hset(
      f"room:{room_id}:geometries",
      "floor",
      json.dumps({"type": Floor.__name__, "data": Floor().model_dump()}),
  )

  3.3. Render Floor in Canvas

  File: app/src/components/Canvas.tsx (line ~72)

  import { Floor } from "./three/Floor";  // ← Add import

  // In geometry mapping (line ~107):
  } else if (config.type === "Floor") {
    return (
      <Floor
        key={name}
        data={config.data}
      />
    );
  }

  ---
  🎯 Key Design Decisions

  1. Static vs Dynamic Data
    - Floor doesn't need frame-based data (unlike Particles)
    - All properties are static configuration
  2. Fog Synchronization
    - Fog automatically syncs with camera.far from settings
    - Uses canvas background color for seamless blend
    - Controlled by fog_density (0-1 range)
  3. Shadow Implementation
    - Using @react-three/drei ContactShadows component
    - Much better performance than Three.js shadow maps
    - Configurable opacity and blur
  4. Grid System
    - Using Three.js GridHelper
    - Dynamic divisions based on size / grid_spacing
    - Customizable colors and opacity
  5. No Base Geometry Features
    - Floor doesn't need selection/hover (not interactive)
    - Overrides position and color as static props
    - No material enum (uses fixed MeshStandardMaterial)

  ---
  📊 Example Configuration

  {
    "type": "Floor",
    "data": {
      "active": true,
      "height": -15.0,
      "size": 200.0,
      "color": "#2a2a2a",
      "grid_spacing": 10.0,
      "grid_color": "#4a4a4a",
      "grid_opacity": 0.6,
      "show_grid": true,
      "show_shadows": true,
      "shadow_opacity": 0.4,
      "shadow_blur": 3.0,
      "fog_enabled": true,
      "fog_density": 0.7
    }
  }

  ---
  ✅ Implementation Order

  1. Backend (Python)
    - Create floor.py with Pydantic model
    - Register in __init__.py
    - Add to /join default geometries
  2. Frontend (TypeScript)
    - Create Floor.tsx component
    - Import in Canvas.tsx
    - Add render case in geometry mapping
