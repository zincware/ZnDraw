// material.tsx
import * as THREE from "three";

/**
 * Material presets for molecular visualization.
 * These are declarative, consistent, and balanced for scientific rendering.
 */
export const MATERIALS = {
  // --- PHYSICAL MATERIALS ---

  "MeshPhysicalMaterial (matt)": {
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

  "MeshPhysicalMaterial (semi-gloss)": {
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

  "MeshPhysicalMaterial (shiny)": {
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

  "MeshPhysicalMaterial (transparent)": {
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

  "MeshPhysicalMaterial (glass)": {
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

  "MeshStandardMaterial (matt)": {
    component: "meshStandardMaterial" as const,
    defaultProps: {
      roughness: 0.8,
      metalness: 0.0,
    },
  },

  "MeshStandardMaterial (metallic)": {
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

  "MeshLambertMaterial (matt)": {
    component: "meshLambertMaterial" as const,
    defaultProps: {},
  },

  "MeshPhongMaterial (classic)": {
    component: "meshPhongMaterial" as const,
    defaultProps: {
      shininess: 30,
      specular: new THREE.Color(0.2, 0.2, 0.2),
    },
  },
} as const;

export const DEFAULT_MATERIAL = "MeshPhysicalMaterial (matt)";
export type MaterialType = keyof typeof MATERIALS;

/**
 * Renders a Three.js material based on a string key.
 * Compatible with Python enum values for consistent cross-language usage.
 */
export function renderMaterial(materialType?: string, opacity?: number, color?: THREE.Color | string, side?: THREE.Side) {
  const type = (materialType || DEFAULT_MATERIAL) as MaterialType;
  const config = MATERIALS[type] || MATERIALS[DEFAULT_MATERIAL];
  const MaterialComponent = config.component;

  const commonProps = {
    color: color || "white",
    side: side || THREE.FrontSide,
    ...config.defaultProps,
    transparent: opacity !== undefined && opacity < 1.0,
    opacity: opacity !== undefined ? opacity : 1.0,
  };

  switch (MaterialComponent) {
    case "meshPhysicalMaterial":
      return <meshPhysicalMaterial {...commonProps} />;
    case "meshStandardMaterial":
      return <meshStandardMaterial {...commonProps} />;
    case "meshBasicMaterial":
      return <meshBasicMaterial {...commonProps} />;
    case "meshToonMaterial":
      return <meshToonMaterial {...commonProps} />;
    case "meshLambertMaterial":
      return <meshLambertMaterial {...commonProps} />;
    case "meshPhongMaterial":
      return <meshPhongMaterial {...commonProps} />;
    default:
      return <meshPhysicalMaterial {...commonProps} />;
  }
}
