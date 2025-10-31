"""Three.js materials module for flexible material configuration.

This module provides Pydantic models for Three.js materials, supporting both
preset strings (e.g., "matt", "shiny") and fully customizable material objects.
"""

from abc import ABC, abstractmethod
from typing import Annotated, Literal

from pydantic import BaseModel, Discriminator, Field


class BaseMaterial(BaseModel, ABC):
    """Base class for all Three.js materials.

    Provides common properties shared across all material types.
    """

    wireframe: bool = Field(default=False, description="Render geometry as wireframe")
    flatShading: bool = Field(
        default=False, description="Use flat shading instead of smooth"
    )
    transparent: bool = Field(default=False, description="Enable transparency")
    polygonOffset: bool = Field(default=False, description="Enable polygon offset")
    polygonOffsetFactor: float = Field(default=0, description="Polygon offset factor")

    @abstractmethod
    def to_three_type(self) -> str:
        """Return the Three.js material type name.

        Returns
        -------
        str
            Three.js material type (e.g., 'meshStandardMaterial')
        """
        pass


class MeshBasicMaterial(BaseMaterial):
    """Simple unlit material."""

    material_type: Literal["MeshBasicMaterial"] = Field(
        default="MeshBasicMaterial", description="Material type discriminator"
    )
    toneMapped: bool = Field(
        default=True, description="Whether material is affected by tone mapping"
    )

    def to_three_type(self) -> str:
        return "meshBasicMaterial"


class MeshStandardMaterial(BaseMaterial):
    """Physically-based rendering (PBR) material.

    Uses metalness/roughness workflow for realistic rendering.
    """

    material_type: Literal["MeshStandardMaterial"] = Field(
        default="MeshStandardMaterial", description="Material type discriminator"
    )
    roughness: float = Field(
        default=1.0, ge=0.0, le=1.0, description="Surface roughness (0=smooth, 1=rough)"
    )
    metalness: float = Field(
        default=0.0,
        ge=0.0,
        le=1.0,
        description="Metallic appearance (0=dielectric, 1=metal)",
    )
    emissive: str = Field(
        default="#000000", description="Emissive (light) color in hex format"
    )
    emissiveIntensity: float = Field(
        default=1.0, ge=0.0, description="Intensity of emissive color"
    )

    def to_three_type(self) -> str:
        return "meshStandardMaterial"


class MeshPhysicalMaterial(MeshStandardMaterial):
    """Advanced physically-based material.

    Extends MeshStandardMaterial with additional physical properties like
    transmission, clearcoat, sheen, and more.
    """

    material_type: Literal["MeshPhysicalMaterial"] = Field(
        default="MeshPhysicalMaterial", description="Material type discriminator"
    )
    transmission: float = Field(
        default=0.0,
        ge=0.0,
        le=1.0,
        description="Light transmission through material (glass effect)",
    )
    thickness: float = Field(
        default=0.0, ge=0.0, description="Thickness for volumetric transmission"
    )
    ior: float = Field(default=1.5, ge=1.0, le=2.333, description="Index of refraction")
    clearcoat: float = Field(
        default=0.0, ge=0.0, le=1.0, description="Clearcoat layer intensity"
    )
    clearcoatRoughness: float = Field(
        default=0.0, ge=0.0, le=1.0, description="Clearcoat layer roughness"
    )
    sheen: float = Field(
        default=0.0, ge=0.0, le=1.0, description="Sheen layer intensity (fabric effect)"
    )
    sheenRoughness: float = Field(
        default=1.0, ge=0.0, le=1.0, description="Sheen layer roughness"
    )
    sheenColor: str = Field(default="#000000", description="Sheen color in hex format")
    specularIntensity: float = Field(
        default=1.0, ge=0.0, description="Specular reflection intensity"
    )
    specularColor: str = Field(
        default="#ffffff", description="Specular reflection color in hex format"
    )
    reflectivity: float = Field(
        default=0.5, ge=0.0, le=1.0, description="Reflectivity intensity"
    )
    envMapIntensity: float = Field(
        default=1.0, ge=0.0, description="Environment map intensity"
    )

    def to_three_type(self) -> str:
        return "meshPhysicalMaterial"


class MeshToonMaterial(BaseMaterial):
    """Non-photorealistic toon/cel-shaded material."""

    material_type: Literal["MeshToonMaterial"] = Field(
        default="MeshToonMaterial", description="Material type discriminator"
    )

    def to_three_type(self) -> str:
        return "meshToonMaterial"


class MeshLambertMaterial(BaseMaterial):
    """Non-shiny diffuse material (Lambertian reflectance)."""

    material_type: Literal["MeshLambertMaterial"] = Field(
        default="MeshLambertMaterial", description="Material type discriminator"
    )
    emissive: str = Field(
        default="#000000", description="Emissive (light) color in hex format"
    )
    emissiveIntensity: float = Field(
        default=1.0, ge=0.0, description="Intensity of emissive color"
    )

    def to_three_type(self) -> str:
        return "meshLambertMaterial"


class MeshPhongMaterial(BaseMaterial):
    """Shiny material with specular highlights (Phong shading)."""

    material_type: Literal["MeshPhongMaterial"] = Field(
        default="MeshPhongMaterial", description="Material type discriminator"
    )
    emissive: str = Field(
        default="#000000", description="Emissive (light) color in hex format"
    )
    emissiveIntensity: float = Field(
        default=1.0, ge=0.0, description="Intensity of emissive color"
    )
    specular: str = Field(
        default="#111111", description="Specular highlight color in hex format"
    )
    shininess: float = Field(
        default=30.0, ge=0.0, description="Shininess/glossiness of specular highlight"
    )

    def to_three_type(self) -> str:
        return "meshPhongMaterial"


# Type union for material properties: preset string or full material object
MaterialType = Annotated[
    MeshBasicMaterial
    | MeshStandardMaterial
    | MeshPhysicalMaterial
    | MeshToonMaterial
    | MeshLambertMaterial
    | MeshPhongMaterial,
    Discriminator("material_type"),
]


# Preset material instances
# These match the TypeScript definitions in app/src/components/three/materials.tsx
PRESETS: dict[str, MaterialType] = {
    # MeshPhysicalMaterial presets
    "MeshPhysicalMaterial_matt": MeshPhysicalMaterial(
        roughness=0.9,
        reflectivity=0.1,
        clearcoat=0.0,
        metalness=0.0,
        transmission=0.0,
        ior=1.45,
    ),
    "MeshPhysicalMaterial_semi-gloss": MeshPhysicalMaterial(
        roughness=0.5,
        reflectivity=0.4,
        clearcoat=0.2,
        metalness=0.0,
        transmission=0.0,
        ior=1.45,
    ),
    "MeshPhysicalMaterial_shiny": MeshPhysicalMaterial(
        roughness=0.2,
        reflectivity=0.6,
        clearcoat=0.4,
        clearcoatRoughness=0.1,
        metalness=0.0,
        transmission=0.0,
        ior=1.45,
    ),
    "MeshPhysicalMaterial_transparent": MeshPhysicalMaterial(
        roughness=0.6,
        reflectivity=0.3,
        clearcoat=0.1,
        transmission=0.6,
        transparent=True,
        thickness=0.2,
        ior=1.2,
    ),
    "MeshPhysicalMaterial_glass": MeshPhysicalMaterial(
        roughness=0.05,
        reflectivity=0.7,
        clearcoat=0.5,
        clearcoatRoughness=0.1,
        transmission=1.0,
        transparent=True,
        thickness=0.5,
        ior=1.5,
        envMapIntensity=1.0,
    ),
    # MeshStandardMaterial presets
    "MeshStandardMaterial_matt": MeshStandardMaterial(
        roughness=0.8,
        metalness=0.0,
    ),
    "MeshStandardMaterial_metallic": MeshStandardMaterial(
        roughness=0.3,
        metalness=1.0,
    ),
    # Simple materials
    "MeshBasicMaterial": MeshBasicMaterial(
        toneMapped=False,
    ),
    "MeshToonMaterial": MeshToonMaterial(),
    "MeshLambertMaterial": MeshLambertMaterial(),
    "MeshPhongMaterial": MeshPhongMaterial(
        shininess=30,
        specular="#333333",  # THREE.Color(0.2, 0.2, 0.2) ≈ #333333
    ),
}


def get_preset(name: str) -> MaterialType | None:
    """Get a material preset by name.

    Parameters
    ----------
    name
        Preset name (e.g., 'MeshPhysicalMaterial_matt')

    Returns
    -------
    MaterialType | None
        The material instance, or None if preset doesn't exist
    """
    return PRESETS.get(name)


def get_all_preset_names() -> list[str]:
    """Get list of all available preset names.

    Returns
    -------
    list[str]
        All available preset names
    """
    return list(PRESETS.keys())


# Type alias for material property: preset string or material object
MaterialProp = (
    Literal[
        "MeshPhysicalMaterial_matt",
        "MeshPhysicalMaterial_semi-gloss",
        "MeshPhysicalMaterial_shiny",
        "MeshPhysicalMaterial_transparent",
        "MeshPhysicalMaterial_glass",
        "MeshStandardMaterial_matt",
        "MeshStandardMaterial_metallic",
        "MeshBasicMaterial",
        "MeshToonMaterial",
        "MeshLambertMaterial",
        "MeshPhongMaterial",
    ]
    | MaterialType
)
