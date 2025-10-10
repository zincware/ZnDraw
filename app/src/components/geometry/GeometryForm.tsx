import { useEffect, useState, useMemo, useRef, useCallback } from "react";
import { debounce } from "lodash";
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
  Chip,
  FormControlLabel,
  Switch,
} from "@mui/material";
import ArrowBackIcon from "@mui/icons-material/ArrowBack";
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
  const clientId = null; // we don't want to skip when saving the current clientId, so undefined
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
  const [activeState, setActiveState] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [isSaving, setIsSaving] = useState(false);
  const isInitializedRef = useRef(false);

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
    isInitializedRef.current = false;
    if (mode === "edit" && geometryData) {
      setKeyInput(geometryData.key);
      setSelectedType(geometryData.geometry.type);
      setFormData(geometryData.geometry.data || {});
      setActiveState(geometryData.geometry.data?.active !== false); // Default to true if undefined
    } else if (mode === "create") {
      setKeyInput("");
      setFormData({});
      setActiveState(true); // Default to true for new geometries
    }
    // Mark as initialized after a short delay to avoid saving initial data
    setTimeout(() => {
      isInitializedRef.current = true;
    }, 500);
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

  const saveGeometry = useCallback(() => {
    if (!roomId || !keyInput.trim() || !selectedType) {
      return;
    }

    setError(null);
    setIsSaving(true);

    createGeometry(
      {
        roomId,
        clientId,
        key: keyInput.trim(),
        geometryType: selectedType,
        geometryData: {
          ...formData,
          active: activeState,
        },
      },
      {
        onSuccess: () => {
          setIsSaving(false);
        },
        onError: (error: any) => {
          setIsSaving(false);
          setError(error.message || "Failed to save geometry");
        },
      }
    );
  }, [roomId, clientId, keyInput, selectedType, formData, activeState, createGeometry]);

  // Create debounced version of save
  const debouncedSave = useMemo(
    () => debounce(saveGeometry, 250),
    [saveGeometry]
  );

  // Auto-save effect
  useEffect(() => {
    // Skip if not initialized or in create mode without key
    if (!isInitializedRef.current || (mode === "create" && !keyInput.trim())) {
      return;
    }

    debouncedSave();

    // Cleanup
    return () => {
      debouncedSave.cancel();
    };
  }, [formData, keyInput, selectedType, mode, debouncedSave]);

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
        <Box sx={{ display: "flex", alignItems: "center", justifyContent: "space-between", mb: 2 }}>
          <Button
            startIcon={<ArrowBackIcon />}
            onClick={handleCancel}
          >
            Back to List
          </Button>
          {isSaving && (
            <Chip
              label="Saving..."
              size="small"
              color="primary"
              sx={{ ml: 2 }}
            />
          )}
        </Box>
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

        <FormControlLabel
          control={
            <Switch
              checked={activeState}
              onChange={(e) => setActiveState(e.target.checked)}
            />
          }
          label="Active (visible in scene)"
          sx={{ mb: 2 }}
        />

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

        <Box sx={{ mt: 3 }}>
          <Button variant="outlined" onClick={handleCancel} fullWidth>
            Close
          </Button>
        </Box>
      </Box>
    </Box>
  );
};

export default GeometryForm;
