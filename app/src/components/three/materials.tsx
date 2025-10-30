// material.tsx
import * as THREE from "three";

/**
 * TypeScript type definitions for Three.js materials.
 * These match the Pydantic models in src/zndraw/materials.py
 */

// Base material properties (common to all materials)
interface BaseMaterialProps {
  wireframe?: boolean;
  flatShading?: boolean;
  transparent?: boolean;
  polygonOffset?: boolean;
  polygonOffsetFactor?: number;
}

export interface MeshBasicMaterialData extends BaseMaterialProps {
  material_type: "MeshBasicMaterial";
  toneMapped?: boolean;
}

export interface MeshStandardMaterialData extends BaseMaterialProps {
  material_type: "MeshStandardMaterial";
  roughness?: number;
  metalness?: number;
  emissive?: string;
  emissiveIntensity?: number;
}

export interface MeshPhysicalMaterialData extends BaseMaterialProps {
  material_type: "MeshPhysicalMaterial";
  roughness?: number;
  metalness?: number;
  emissive?: string;
  emissiveIntensity?: number;
  transmission?: number;
  thickness?: number;
  ior?: number;
  reflectivity?: number;
  clearcoat?: number;
  clearcoatRoughness?: number;
  sheen?: number;
  sheenRoughness?: number;
  sheenColor?: string;
  specularIntensity?: number;
  specularColor?: string;
  envMapIntensity?: number;
}

export interface MeshToonMaterialData extends BaseMaterialProps {
  material_type: "MeshToonMaterial";
}

export interface MeshLambertMaterialData extends BaseMaterialProps {
  material_type: "MeshLambertMaterial";
  emissive?: string;
  emissiveIntensity?: number;
}

export interface MeshPhongMaterialData extends BaseMaterialProps {
  material_type: "MeshPhongMaterial";
  emissive?: string;
  emissiveIntensity?: number;
  specular?: string;
  shininess?: number;
}

export type MaterialData =
  | MeshBasicMaterialData
  | MeshStandardMaterialData
  | MeshPhysicalMaterialData
  | MeshToonMaterialData
  | MeshLambertMaterialData
  | MeshPhongMaterialData;

export type MaterialProp = string | MaterialData;

/**
 * Material presets for backward compatibility.
 * Note: These are now defined in Python (src/zndraw/materials.py) as the source of truth.
 * This mapping is kept for fallback when string presets are used.
 */
export const MATERIALS = {
  // --- PHYSICAL MATERIALS ---

  "MeshPhysicalMaterial_matt": {
    component: "meshPhysicalMaterial" as const,
    defaultProps: {
      roughness: 0.9,           // very diffuse
      reflectivity: 0.1,        // minimal reflections
      clearcoat: 0.0,
      metalness: 0.0,
      transmission: 0.0,        // fully opaque
      ior: 1.45,
    },
  },

  "MeshPhysicalMaterial_semi-gloss": {
    component: "meshPhysicalMaterial" as const,
    defaultProps: {
      roughness: 0.5,
      reflectivity: 0.4,
      clearcoat: 0.2,
      metalness: 0.0,
      transmission: 0.0,
      ior: 1.45,
    },
  },

  "MeshPhysicalMaterial_shiny": {
    component: "meshPhysicalMaterial" as const,
    defaultProps: {
      roughness: 0.2,
      reflectivity: 0.6,
      clearcoat: 0.4,
      clearcoatRoughness: 0.1,
      metalness: 0.0,
      transmission: 0.0,
      ior: 1.45,
    },
  },

  "MeshPhysicalMaterial_transparent": {
    component: "meshPhysicalMaterial" as const,
    defaultProps: {
      roughness: 0.6,
      reflectivity: 0.3,
      clearcoat: 0.1,
      transmission: 0.6,   // light passes through partially
      opacity: 0.6,        // visible transparency strength
      transparent: true,   // required for opacity/transmission
      thickness: 0.2,
      ior: 1.2,
      depthWrite: false,   // avoid depth conflicts when layered
    },
  },

  "MeshPhysicalMaterial_glass": {
    component: "meshPhysicalMaterial" as const,
    defaultProps: {
      roughness: 0.05,
      reflectivity: 0.7,
      clearcoat: 0.5,
      clearcoatRoughness: 0.1,
      transmission: 1.0,   // fully light-transmitting
      opacity: 0.25,       // partial visibility; tweak as desired
      transparent: true,   // essential for opacity
      thickness: 0.5,
      ior: 1.5,            // realistic glass IOR
      envMapIntensity: 1.0,
      depthWrite: false,
    },
  },

  // --- STANDARD MATERIALS ---

  "MeshStandardMaterial_matt": {
    component: "meshStandardMaterial" as const,
    defaultProps: {
      roughness: 0.8,
      metalness: 0.0,
    },
  },

  "MeshStandardMaterial_metallic": {
    component: "meshStandardMaterial" as const,
    defaultProps: {
      roughness: 0.3,
      metalness: 1.0,
    },
  },

  // --- SIMPLE / STYLISED MATERIALS ---

  "MeshBasicMaterial": {
    component: "meshBasicMaterial" as const,
    defaultProps: {
      toneMapped: false, // unaffected by lighting — for overlays/debug
    },
  },

  "MeshToonMaterial": {
    component: "meshToonMaterial" as const,
    defaultProps: {
      gradientMap: undefined,
    },
  },

  "MeshLambertMaterial": {
    component: "meshLambertMaterial" as const,
    defaultProps: {},
  },

  "MeshPhongMaterial": {
    component: "meshPhongMaterial" as const,
    defaultProps: {
      shininess: 30,
      specular: new THREE.Color(0.2, 0.2, 0.2),
    },
  },
} as const;

export const DEFAULT_MATERIAL = "MeshPhysicalMaterial_matt";
export type MaterialType = keyof typeof MATERIALS;

/**
 * Type guard to check if material is a string preset
 */
function isStringMaterial(material: MaterialProp | undefined): material is string {
  return typeof material === "string";
}

/**
 * Type guard to check if material is a material object
 */
function isMaterialObject(material: MaterialProp | undefined): material is MaterialData {
  return typeof material === "object" && material !== null && "material_type" in material;
}

/**
 * Renders a Three.js material based on either a string preset or material object.
 * Supports both legacy string presets and new material objects from Python.
 *
 * @param material - Material preset string or material object
 * @param opacity - Geometry opacity (overrides material opacity)
 * @param color - Geometry color (overrides material color)
 * @param side - Rendering side (auto-determined if not specified)
 */
export function renderMaterial(
  material?: MaterialProp,
  opacity?: number,
  color?: THREE.Color | string,
  side?: THREE.Side
) {
  // Common properties applied to all materials
  const commonProps = {
    color: color || "white",
    side: side || (opacity !== undefined && opacity < 1.0 ? THREE.DoubleSide : THREE.FrontSide),
    transparent: opacity !== undefined && opacity < 1.0,
    opacity: opacity !== undefined ? opacity : 1.0,
  };

  // Handle material object
  if (isMaterialObject(material)) {
    // Extract properties from material object
    const materialProps = {
      ...commonProps,
      wireframe: material.wireframe,
      flatShading: material.flatShading,
      transparent: material.transparent || commonProps.transparent,
      polygonOffset: material.polygonOffset,
      polygonOffsetFactor: material.polygonOffsetFactor,
    };

    switch (material.material_type) {
      case "MeshBasicMaterial":
        return (
          <meshBasicMaterial
            {...materialProps}
            toneMapped={material.toneMapped}
          />
        );

      case "MeshStandardMaterial":
        return (
          <meshStandardMaterial
            {...materialProps}
            roughness={material.roughness}
            metalness={material.metalness}
            emissive={material.emissive}
            emissiveIntensity={material.emissiveIntensity}
          />
        );

      case "MeshPhysicalMaterial":
        return (
          <meshPhysicalMaterial
            {...materialProps}
            roughness={material.roughness}
            metalness={material.metalness}
            emissive={material.emissive}
            emissiveIntensity={material.emissiveIntensity}
            transmission={material.transmission}
            thickness={material.thickness}
            ior={material.ior}
            reflectivity={material.reflectivity}
            clearcoat={material.clearcoat}
            clearcoatRoughness={material.clearcoatRoughness}
            sheen={material.sheen}
            sheenRoughness={material.sheenRoughness}
            sheenColor={material.sheenColor}
            specularIntensity={material.specularIntensity}
            specularColor={material.specularColor}
            envMapIntensity={material.envMapIntensity}
          />
        );

      case "MeshToonMaterial":
        return <meshToonMaterial {...materialProps} />;

      case "MeshLambertMaterial":
        return (
          <meshLambertMaterial
            {...materialProps}
            emissive={material.emissive}
            emissiveIntensity={material.emissiveIntensity}
          />
        );

      case "MeshPhongMaterial":
        return (
          <meshPhongMaterial
            {...materialProps}
            emissive={material.emissive}
            emissiveIntensity={material.emissiveIntensity}
            specular={material.specular}
            shininess={material.shininess}
          />
        );
    }
  }

  // Handle string preset (backward compatibility)
  if (isStringMaterial(material)) {
    const type = material as MaterialType;
    const config = MATERIALS[type] || MATERIALS[DEFAULT_MATERIAL];
    const MaterialComponent = config.component;

    const presetProps = {
      ...commonProps,
      ...config.defaultProps,
    };

    switch (MaterialComponent) {
      case "meshPhysicalMaterial":
        return <meshPhysicalMaterial {...presetProps} />;
      case "meshStandardMaterial":
        return <meshStandardMaterial {...presetProps} />;
      case "meshBasicMaterial":
        return <meshBasicMaterial {...presetProps} />;
      case "meshToonMaterial":
        return <meshToonMaterial {...presetProps} />;
      case "meshLambertMaterial":
        return <meshLambertMaterial {...presetProps} />;
      case "meshPhongMaterial":
        return <meshPhongMaterial {...presetProps} />;
      default:
        return <meshPhysicalMaterial {...presetProps} />;
    }
  }

  // Default fallback
  return <meshPhysicalMaterial {...commonProps} roughness={0.9} reflectivity={0.1} />;
}
