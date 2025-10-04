import * as THREE from "three";

/**
 * Available material types for geometries.
 * Each material has a component name and default props.
 */
export const MATERIALS = {
  MeshPhysicalMaterial: {
    component: "meshPhysicalMaterial" as const,
    defaultProps: {
      roughness: 0.3,
      reflectivity: 0.4,
      clearcoat: 0.1,
    },
  },
  MeshStandardMaterial: {
    component: "meshStandardMaterial" as const,
    defaultProps: {
      roughness: 0.5,
      metalness: 0.0,
    },
  },
  MeshBasicMaterial: {
    component: "meshBasicMaterial" as const,
    defaultProps: {},
  },
  MeshToonMaterial: {
    component: "meshToonMaterial" as const,
    defaultProps: {},
  },
} as const;

export const DEFAULT_MATERIAL = "MeshPhysicalMaterial";

export type MaterialType = keyof typeof MATERIALS;

/**
 * Renders the appropriate Three.js material based on material type string.
 *
 * @param materialType - The material type name (e.g., "MeshPhysicalMaterial")
 * @returns JSX element for the material component
 */
export function renderMaterial(materialType?: string) {
  const type = (materialType || DEFAULT_MATERIAL) as MaterialType;
  const config = MATERIALS[type] || MATERIALS[DEFAULT_MATERIAL];
  const MaterialComponent = config.component;

  const commonProps = {
    color: "white",
    side: THREE.FrontSide,
    ...config.defaultProps,
  };

  // TypeScript needs explicit handling for each component type
  switch (MaterialComponent) {
    case "meshPhysicalMaterial":
      return <meshPhysicalMaterial {...commonProps} />;
    case "meshStandardMaterial":
      return <meshStandardMaterial {...commonProps} />;
    case "meshBasicMaterial":
      return <meshBasicMaterial {...commonProps} />;
    case "meshToonMaterial":
      return <meshToonMaterial {...commonProps} />;
    default:
      return <meshPhysicalMaterial {...commonProps} />;
  }
}
