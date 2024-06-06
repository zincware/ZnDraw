import enum

from pydantic import BaseModel, Field


class Material(str, enum.Enum):
    MeshBasicMaterial = "MeshBasicMaterial"
    MeshLambertMaterial = "MeshLambertMaterial"
    MeshMatcapMaterial = "MeshMatcapMaterial"
    MeshPhongMaterial = "MeshPhongMaterial"
    MeshPhysicalMaterial = "MeshPhysicalMaterial"
    MeshStandardMaterial = "MeshStandardMaterial"
    MeshToonMaterial = "MeshToonMaterial"


# create a class for the material, resolution, etc.
class Scene(BaseModel):
    fps: int = Field(30, ge=1, le=120, description="Maxium frames per second")
    # material: Material = Field(Material.MeshPhongMaterial, description="Material")
    # resolution: int = Field(10, ge=1, le=50, description="Resolution")
    # particle_size: float = Field(1.0, ge=0.1, le=5, description="Particle Size")
    # bonds_size: float = Field(1.0, ge=0.1, le=5, description="Bonds Size")
    # wireframe: bool = Field(False, description="Wireframe")
    loop: bool = Field(
        False,
        alias="Animation Loop",
        description="Automatically restart animation when finished.",
    )
    simulation_box: bool = Field(
        False,
        description="Show the simulation box.",
    )
    # bonds: bool = Field(
    #     True,
    #     description="Show bonds.",
    # )
    # line_label: bool = Field(
    #     True,
    #     description="Show the length of the line.",
    # )
    # label_offset: int = Field(
    #     0,
    #     ge=-7,
    #     le=7,
    #     description="Move the label to the left or right (keypress i).",
    # )

    @classmethod
    def updated_schema(cls) -> dict:
        schema = cls.model_json_schema()

        # schema["properties"]["wireframe"]["format"] = "checkbox"
        schema["properties"]["Animation Loop"]["format"] = "checkbox"
        schema["properties"]["simulation_box"]["format"] = "checkbox"
        # schema["properties"]["resolution"]["format"] = "range"
        # schema["properties"]["label_offset"]["format"] = "range"
        # schema["properties"]["particle_size"]["format"] = "range"
        schema["properties"]["fps"]["format"] = "range"
        # schema["properties"]["particle_size"]["step"] = 0.1
        # schema["properties"]["bonds_size"]["format"] = "range"
        # schema["properties"]["bonds_size"]["step"] = 0.1
        # schema["properties"]["bonds"]["format"] = "checkbox"
        # schema["properties"]["line_label"]["format"] = "checkbox"

        return schema
