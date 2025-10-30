import { useState } from "react";
import { withJsonFormsControlProps } from "@jsonforms/react";
import { rankWith, schemaMatches, and, uiTypeIs, ControlProps } from "@jsonforms/core";
import {
  Box,
  Button,
  FormControl,
  FormLabel,
  MenuItem,
  Select,
  TextField,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Typography,
  Checkbox,
  FormControlLabel,
  Stack,
  Chip,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import EditIcon from "@mui/icons-material/Edit";
import RestartAltIcon from "@mui/icons-material/RestartAlt";
import { MaterialData, MaterialProp } from "../three/materials";

/**
 * Material presets matching Python definitions
 */
const MATERIAL_PRESETS = [
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
];

/**
 * Material types for the type selector
 */
const MATERIAL_TYPES = [
  { value: "MeshBasicMaterial", label: "Basic (Unlit)" },
  { value: "MeshStandardMaterial", label: "Standard (PBR)" },
  { value: "MeshPhysicalMaterial", label: "Physical (Advanced PBR)" },
  { value: "MeshToonMaterial", label: "Toon (Cel-Shaded)" },
  { value: "MeshLambertMaterial", label: "Lambert (Diffuse)" },
  { value: "MeshPhongMaterial", label: "Phong (Specular)" },
];

/**
 * Type guard to check if value is a string preset
 */
function isStringMaterial(value: MaterialProp | undefined): value is string {
  return typeof value === "string";
}

/**
 * Type guard to check if value is a material object
 */
function isMaterialObject(value: MaterialProp | undefined): value is MaterialData {
  return typeof value === "object" && value !== null && "material_type" in value;
}

/**
 * Create a default material object for a given type
 */
function createDefaultMaterial(materialType: MaterialData["material_type"]): MaterialData {
  const base = {
    wireframe: false,
    flatShading: false,
    transparent: false,
    polygonOffset: false,
    polygonOffsetFactor: 0,
  };

  switch (materialType) {
    case "MeshBasicMaterial":
      return { ...base, material_type: "MeshBasicMaterial", toneMapped: true };

    case "MeshStandardMaterial":
      return {
        ...base,
        material_type: "MeshStandardMaterial",
        roughness: 1.0,
        metalness: 0.0,
        emissive: "#000000",
        emissiveIntensity: 1.0,
      };

    case "MeshPhysicalMaterial":
      return {
        ...base,
        material_type: "MeshPhysicalMaterial",
        roughness: 1.0,
        metalness: 0.0,
        emissive: "#000000",
        emissiveIntensity: 1.0,
        transmission: 0.0,
        thickness: 0.0,
        ior: 1.5,
        reflectivity: 0.5,
        clearcoat: 0.0,
        clearcoatRoughness: 0.0,
        sheen: 0.0,
        sheenRoughness: 1.0,
        sheenColor: "#000000",
        specularIntensity: 1.0,
        specularColor: "#ffffff",
        envMapIntensity: 1.0,
      };

    case "MeshToonMaterial":
      return { ...base, material_type: "MeshToonMaterial" };

    case "MeshLambertMaterial":
      return {
        ...base,
        material_type: "MeshLambertMaterial",
        emissive: "#000000",
        emissiveIntensity: 1.0,
      };

    case "MeshPhongMaterial":
      return {
        ...base,
        material_type: "MeshPhongMaterial",
        emissive: "#000000",
        emissiveIntensity: 1.0,
        specular: "#111111",
        shininess: 30.0,
      };
  }
}

/**
 * MaterialEditor Component
 *
 * Provides a hybrid UI for material editing:
 * - Preset dropdown + Edit button (when value is string)
 * - Full property editor (when value is object)
 * - Reset to Preset button (to convert object back to string)
 */
const MaterialEditor = ({
  data,
  handleChange,
  path,
  label,
  errors,
}: ControlProps) => {
  const [expanded, setExpanded] = useState(false);

  const isPreset = isStringMaterial(data);
  const isObject = isMaterialObject(data);

  /**
   * Convert preset string to material object for editing
   */
  const handleConvertToObject = () => {
    if (isStringMaterial(data)) {
      // Create a default material object (user can customize from here)
      const defaultMaterial = createDefaultMaterial("MeshPhysicalMaterial");
      handleChange(path, defaultMaterial);
      setExpanded(true);
    }
  };

  /**
   * Reset to a preset string
   */
  const handleResetToPreset = () => {
    handleChange(path, "MeshPhysicalMaterial_matt");
    setExpanded(false);
  };

  /**
   * Update material object property
   */
  const handleUpdateProperty = (property: string, value: any) => {
    if (isMaterialObject(data)) {
      handleChange(path, { ...data, [property]: value });
    }
  };

  /**
   * Change material type (creates new object with default values for that type)
   */
  const handleChangeMaterialType = (newType: MaterialData["material_type"]) => {
    const newMaterial = createDefaultMaterial(newType);
    // Preserve common properties if they exist
    if (isMaterialObject(data)) {
      newMaterial.wireframe = data.wireframe;
      newMaterial.flatShading = data.flatShading;
      newMaterial.transparent = data.transparent;
      newMaterial.polygonOffset = data.polygonOffset;
      newMaterial.polygonOffsetFactor = data.polygonOffsetFactor;
    }
    handleChange(path, newMaterial);
  };

  return (
    <Box sx={{ mb: 2 }}>
      <FormControl fullWidth error={!!errors}>
        <FormLabel sx={{ mb: 1 }}>
          {label || "Material"}
          {isObject && (
            <Chip
              label={data.material_type}
              size="small"
              sx={{ ml: 1 }}
              color="primary"
              variant="outlined"
            />
          )}
        </FormLabel>

        {/* Preset Dropdown (when value is string) */}
        {isPreset && (
          <Stack direction="row" spacing={1} alignItems="center">
            <Select
              value={data || "MeshPhysicalMaterial_matt"}
              onChange={(e) => handleChange(path, e.target.value)}
              fullWidth
              size="small"
            >
              {MATERIAL_PRESETS.map((preset) => (
                <MenuItem key={preset} value={preset}>
                  {preset.replace(/_/g, " ")}
                </MenuItem>
              ))}
            </Select>
            <Button
              variant="outlined"
              startIcon={<EditIcon />}
              onClick={handleConvertToObject}
              size="small"
            >
              Edit
            </Button>
          </Stack>
        )}

        {/* Material Object Editor (when value is object) */}
        {isObject && (
          <Box>
            {/* Material Type Selector */}
            <Box sx={{ mb: 2 }}>
              <FormControl fullWidth size="small">
                <FormLabel sx={{ mb: 0.5, fontSize: "0.875rem" }}>
                  Material Type
                </FormLabel>
                <Select
                  value={data.material_type}
                  onChange={(e) =>
                    handleChangeMaterialType(
                      e.target.value as MaterialData["material_type"]
                    )
                  }
                  size="small"
                >
                  {MATERIAL_TYPES.map((type) => (
                    <MenuItem key={type.value} value={type.value}>
                      {type.label}
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
            </Box>

            {/* Common Properties Accordion */}
            <Accordion expanded={expanded} onChange={() => setExpanded(!expanded)}>
              <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                <Typography variant="body2">Material Properties</Typography>
              </AccordionSummary>
              <AccordionDetails>
                <Stack spacing={2}>
                  {/* Common Properties */}
                  <Box>
                    <Typography variant="subtitle2" sx={{ mb: 1 }}>
                      Common
                    </Typography>
                    <Stack spacing={1}>
                      <FormControlLabel
                        control={
                          <Checkbox
                            checked={data.wireframe || false}
                            onChange={(e) =>
                              handleUpdateProperty("wireframe", e.target.checked)
                            }
                            size="small"
                          />
                        }
                        label="Wireframe"
                      />
                      <FormControlLabel
                        control={
                          <Checkbox
                            checked={data.flatShading || false}
                            onChange={(e) =>
                              handleUpdateProperty("flatShading", e.target.checked)
                            }
                            size="small"
                          />
                        }
                        label="Flat Shading"
                      />
                      <FormControlLabel
                        control={
                          <Checkbox
                            checked={data.transparent || false}
                            onChange={(e) =>
                              handleUpdateProperty("transparent", e.target.checked)
                            }
                            size="small"
                          />
                        }
                        label="Transparent"
                      />
                    </Stack>
                  </Box>

                  {/* Material-Specific Properties */}
                  {data.material_type === "MeshStandardMaterial" && (
                    <Box>
                      <Typography variant="subtitle2" sx={{ mb: 1 }}>
                        Standard Material
                      </Typography>
                      <Stack spacing={1.5}>
                        <TextField
                          label="Roughness"
                          type="number"
                          value={(data as any).roughness ?? 1.0}
                          onChange={(e) =>
                            handleUpdateProperty("roughness", parseFloat(e.target.value))
                          }
                          inputProps={{ min: 0, max: 1, step: 0.1 }}
                          size="small"
                          fullWidth
                        />
                        <TextField
                          label="Metalness"
                          type="number"
                          value={(data as any).metalness ?? 0.0}
                          onChange={(e) =>
                            handleUpdateProperty("metalness", parseFloat(e.target.value))
                          }
                          inputProps={{ min: 0, max: 1, step: 0.1 }}
                          size="small"
                          fullWidth
                        />
                      </Stack>
                    </Box>
                  )}

                  {data.material_type === "MeshPhysicalMaterial" && (
                    <Box>
                      <Typography variant="subtitle2" sx={{ mb: 1 }}>
                        Physical Material
                      </Typography>
                      <Stack spacing={1.5}>
                        <TextField
                          label="Roughness"
                          type="number"
                          value={(data as any).roughness ?? 1.0}
                          onChange={(e) =>
                            handleUpdateProperty("roughness", parseFloat(e.target.value))
                          }
                          inputProps={{ min: 0, max: 1, step: 0.1 }}
                          size="small"
                          fullWidth
                        />
                        <TextField
                          label="Metalness"
                          type="number"
                          value={(data as any).metalness ?? 0.0}
                          onChange={(e) =>
                            handleUpdateProperty("metalness", parseFloat(e.target.value))
                          }
                          inputProps={{ min: 0, max: 1, step: 0.1 }}
                          size="small"
                          fullWidth
                        />
                        <TextField
                          label="Transmission"
                          type="number"
                          value={(data as any).transmission ?? 0.0}
                          onChange={(e) =>
                            handleUpdateProperty("transmission", parseFloat(e.target.value))
                          }
                          inputProps={{ min: 0, max: 1, step: 0.1 }}
                          size="small"
                          fullWidth
                        />
                        <TextField
                          label="IOR (Index of Refraction)"
                          type="number"
                          value={(data as any).ior ?? 1.5}
                          onChange={(e) =>
                            handleUpdateProperty("ior", parseFloat(e.target.value))
                          }
                          inputProps={{ min: 1, max: 2.333, step: 0.1 }}
                          size="small"
                          fullWidth
                        />
                        <TextField
                          label="Clearcoat"
                          type="number"
                          value={(data as any).clearcoat ?? 0.0}
                          onChange={(e) =>
                            handleUpdateProperty("clearcoat", parseFloat(e.target.value))
                          }
                          inputProps={{ min: 0, max: 1, step: 0.1 }}
                          size="small"
                          fullWidth
                        />
                      </Stack>
                    </Box>
                  )}

                  {data.material_type === "MeshPhongMaterial" && (
                    <Box>
                      <Typography variant="subtitle2" sx={{ mb: 1 }}>
                        Phong Material
                      </Typography>
                      <Stack spacing={1.5}>
                        <TextField
                          label="Shininess"
                          type="number"
                          value={(data as any).shininess ?? 30.0}
                          onChange={(e) =>
                            handleUpdateProperty("shininess", parseFloat(e.target.value))
                          }
                          inputProps={{ min: 0, step: 1 }}
                          size="small"
                          fullWidth
                        />
                        <TextField
                          label="Specular Color"
                          type="text"
                          value={(data as any).specular ?? "#111111"}
                          onChange={(e) =>
                            handleUpdateProperty("specular", e.target.value)
                          }
                          size="small"
                          fullWidth
                        />
                      </Stack>
                    </Box>
                  )}
                </Stack>
              </AccordionDetails>
            </Accordion>

            {/* Reset to Preset Button */}
            <Box sx={{ mt: 2 }}>
              <Button
                variant="outlined"
                startIcon={<RestartAltIcon />}
                onClick={handleResetToPreset}
                size="small"
                fullWidth
              >
                Reset to Preset
              </Button>
            </Box>
          </Box>
        )}

        {errors && (
          <Typography variant="caption" color="error" sx={{ mt: 0.5 }}>
            {errors}
          </Typography>
        )}
      </FormControl>
    </Box>
  );
};

/**
 * Tester function for MaterialEditor
 * Matches fields with x-custom-type="three-material"
 * Priority 10 to override default renderers
 */
export const materialEditorTester = rankWith(
  10,
  and(
    schemaMatches((schema) => (schema as any)["x-custom-type"] === "three-material"),
    uiTypeIs("Control")
  )
);

export default withJsonFormsControlProps(MaterialEditor);
