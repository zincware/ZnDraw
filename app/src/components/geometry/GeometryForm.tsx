import { useEffect, useState, useMemo } from "react";
import {
  Box,
  Button,
  TextField,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Typography,
  CircularProgress,
  Alert,
  SelectChangeEvent,
} from "@mui/material";
import ArrowBackIcon from "@mui/icons-material/ArrowBack";
import SaveIcon from "@mui/icons-material/Save";
import { JsonForms } from "@jsonforms/react";
import { materialCells } from "@jsonforms/material-renderers";
import { useGeometryStore } from "../../stores/geometryStore";
import {
  useGeometry,
  useGeometrySchemas,
  useCreateGeometry,
} from "../../hooks/useGeometries";
import { useAppStore } from "../../store";
import { useFrameMetadata } from "../../hooks/useSchemas";
import { customRenderers, injectDynamicEnums } from "../../utils/jsonforms";

const GeometryForm = () => {
  const { roomId } = useAppStore();
  const clientId = undefined; // we don't want to skip when saving the current clientId, so undefined
  const {
    mode,
    selectedKey,
    selectedType,
    formData,
    setMode,
    setSelectedType,
    setFormData,
    resetForm,
  } = useGeometryStore();

  const [keyInput, setKeyInput] = useState("");
  const [error, setError] = useState<string | null>(null);

  // Fetch geometry schemas
  const {
    data: schemasData,
    isLoading: isLoadingSchemas,
    isError: isSchemasError,
  } = useGeometrySchemas(roomId);

  // Fetch frame metadata for dynamic enums
  const { data: metadata, isLoading: isLoadingMetadata } =
    useFrameMetadata(roomId);

  // Fetch existing geometry if in edit mode
  const {
    data: geometryData,
    isLoading: isLoadingGeometry,
    isError: isGeometryError,
  } = useGeometry(roomId, mode === "edit" ? selectedKey : null);

  // Create geometry mutation
  const { mutate: createGeometry, isPending: isCreating } =
    useCreateGeometry();

  // Extract schemas from the response
  const schemas = useMemo(() => {
    if (!schemasData?.schemas) return {};
    return schemasData.schemas;
  }, [schemasData]);

  const schemaOptions = useMemo(() => Object.keys(schemas), [schemas]);

  // Initialize form data for edit mode
  useEffect(() => {
    if (mode === "edit" && geometryData) {
      setKeyInput(geometryData.key);
      setSelectedType(geometryData.geometry.type);
      setFormData(geometryData.geometry.data || {});
    } else if (mode === "create") {
      setKeyInput("");
      setFormData({});
    }
  }, [mode, geometryData, setFormData, setSelectedType]);

  const handleTypeChange = (event: SelectChangeEvent<string>) => {
    const type = event.target.value;
    setSelectedType(type);
    setFormData({}); // Reset form data when type changes
    setError(null);
  };

  const handleFormChange = ({ data }: { data: any }) => {
    setFormData(data ?? {});
  };

  const handleSave = () => {
    if (!roomId) {
      setError("Room ID is missing");
      return;
    }

    if (!keyInput.trim()) {
      setError("Key is required");
      return;
    }

    if (!selectedType) {
      setError("Type is required");
      return;
    }

    setError(null);

    createGeometry(
      {
        roomId,
        clientId,
        key: keyInput.trim(),
        geometryType: selectedType,
        geometryData: formData,
      },
      {
        onSuccess: () => {
          resetForm();
        },
        onError: (error: any) => {
          setError(error.message || "Failed to save geometry");
        },
      }
    );
  };

  const handleCancel = () => {
    resetForm();
    setError(null);
  };

  // Create dynamic schema with injected metadata
  const currentSchema = useMemo(() => {
    if (!selectedType || !schemas[selectedType]) return null;
    return injectDynamicEnums(schemas[selectedType], metadata);
  }, [selectedType, schemas, metadata]);

  if (isLoadingSchemas || isLoadingMetadata) {
    return (
      <Box
        sx={{
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
          height: "100%",
        }}
      >
        <CircularProgress />
      </Box>
    );
  }

  if (isSchemasError) {
    return (
      <Box sx={{ p: 2 }}>
        <Typography color="error">Failed to load geometry schemas.</Typography>
        <Button onClick={handleCancel} sx={{ mt: 2 }}>
          Back to List
        </Button>
      </Box>
    );
  }

  if (mode === "edit" && (isLoadingGeometry || isGeometryError)) {
    return (
      <Box sx={{ p: 2 }}>
        {isLoadingGeometry && <CircularProgress />}
        {isGeometryError && (
          <>
            <Typography color="error">Failed to load geometry.</Typography>
            <Button onClick={handleCancel} sx={{ mt: 2 }}>
              Back to List
            </Button>
          </>
        )}
      </Box>
    );
  }

  return (
    <Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
      <Box sx={{ p: 2, borderBottom: "1px solid rgba(0, 0, 0, 0.12)" }}>
        <Button
          startIcon={<ArrowBackIcon />}
          onClick={handleCancel}
          sx={{ mb: 2 }}
        >
          Back to List
        </Button>
        <Typography variant="h6">
          {mode === "create" ? "Create Geometry" : "Edit Geometry"}
        </Typography>
      </Box>

      <Box sx={{ p: 2, pb: 12, flexGrow: 1, overflowY: "auto" }}>
        {error && (
          <Alert severity="error" sx={{ mb: 2 }} onClose={() => setError(null)}>
            {error}
          </Alert>
        )}

        <TextField
          fullWidth
          label="Key"
          value={keyInput}
          onChange={(e) => setKeyInput(e.target.value)}
          disabled={mode === "edit"}
          sx={{ mb: 2 }}
          required
        />

        <FormControl fullWidth sx={{ mb: 2 }} required>
          <InputLabel>Type</InputLabel>
          <Select
            value={selectedType || ""}
            label="Type"
            onChange={handleTypeChange}
            disabled={mode === "edit"}
          >
            {schemaOptions.map((type) => (
              <MenuItem key={type} value={type}>
                {type}
              </MenuItem>
            ))}
          </Select>
        </FormControl>

        {currentSchema && (
          <>
            <Typography variant="subtitle2" sx={{ mb: 2 }}>
              Geometry Data
            </Typography>
            <JsonForms
              schema={currentSchema}
              data={formData}
              renderers={customRenderers}
              cells={materialCells}
              onChange={handleFormChange}
            />
          </>
        )}

        <Box sx={{ mt: 3, display: "flex", gap: 2 }}>
          <Button
            variant="contained"
            startIcon={<SaveIcon />}
            onClick={handleSave}
            disabled={isCreating || !selectedType || !keyInput.trim()}
            fullWidth
          >
            {isCreating ? "Saving..." : "Save"}
          </Button>
          <Button variant="outlined" onClick={handleCancel} fullWidth>
            Cancel
          </Button>
        </Box>
      </Box>
    </Box>
  );
};

export default GeometryForm;
