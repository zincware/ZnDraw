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
import { useFrameKeys } from "../../hooks/useSchemas";
import { customRenderers, injectDynamicEnums } from "../../utils/jsonforms";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

const GeometryForm = () => {
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const geometries = useAppStore((state) => state.geometries);
  const geometryDefaults = useAppStore((state) => state.geometryDefaults);
  const userName = null; // we don't want to skip when saving the current userName, so undefined
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
  const [isLocked, setIsLocked] = useState(false); // Track if key/type are locked
  const isInitializedRef = useRef(false);

  // Fetch geometry schemas
  const {
    data: schemasData,
    isLoading: isLoadingSchemas,
    isError: isSchemasError,
  } = useGeometrySchemas(roomId);

  // Fetch frame keys for dynamic enums
  const { data: metadata, isLoading: isLoadingMetadata } =
    useFrameKeys(roomId!);

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
      isInitializedRef.current = false;
      setKeyInput(geometryData.key);
      setSelectedType(geometryData.geometry.type);
      setFormData(geometryData.geometry.data || {});
      setActiveState(geometryData.geometry.data?.active !== false); // Default to true if undefined
      setIsLocked(true); // Lock fields in edit mode

      // Mark as initialized after a short delay to avoid saving initial data
      setTimeout(() => {
        isInitializedRef.current = true;
      }, 500);
    }
  }, [mode, geometryData, setFormData, setSelectedType]);

  // Apply defaults and lock fields when both key and type are present in create mode
  const lockAndCreateGeometry = useCallback(() => {
    if (mode === "create" && keyInput.trim() && selectedType && geometryDefaults && !isLocked) {
      const defaultData = getGeometryWithDefaults({}, selectedType, geometryDefaults);
      setFormData(defaultData);
      setIsLocked(true);
      isInitializedRef.current = true;
    }
  }, [mode, keyInput, selectedType, geometryDefaults, isLocked, setFormData]);

  // Handler for when key or type field is blurred
  const handleFieldBlur = useCallback(() => {
    if (keyInput.trim() && selectedType) {
      lockAndCreateGeometry();
    }
  }, [keyInput, selectedType, lockAndCreateGeometry]);

  const handleTypeChange = (event: SelectChangeEvent<string>) => {
    const type = event.target.value;
    setSelectedType(type);
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
          const errorMsg = error.message || "Failed to save geometry";
          setError(errorMsg);
          // Snackbar is now shown automatically by withAutoLock for lock failures
        },
      }
    );
  }, [roomId, keyInput, selectedType, formData, activeState, createGeometry]);

  // Create debounced version of save
  const debouncedSave = useMemo(
    () => debounce(saveGeometry, 250),
    [saveGeometry]
  );

  // Auto-save effect
  useEffect(() => {
    // Skip if not initialized or in create mode without key/type
    if (!isInitializedRef.current || (mode === "create" && (!keyInput.trim() || !selectedType))) {
      return;
    }

    debouncedSave();

    // Cleanup
    return () => {
      debouncedSave.cancel();
    };
  }, [formData, keyInput, selectedType, activeState, mode, debouncedSave]);

  const handleCancel = () => {
    resetForm();
    setError(null);
    setIsLocked(false);
    setKeyInput("");
    isInitializedRef.current = false;
  };

  // Create dynamic schema with injected metadata and geometries
  const currentSchema = useMemo(() => {
    if (!selectedType || !schemas[selectedType]) return null;
    return injectDynamicEnums(schemas[selectedType], metadata, geometries);
  }, [selectedType, schemas, metadata, geometries]);

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
          onBlur={handleFieldBlur}
          disabled={isLocked}
          sx={{ mb: 2 }}
          required
          helperText={!isLocked && mode === "create" ? "Enter a unique name for this geometry" : ""}
        />

        <FormControl fullWidth sx={{ mb: 2 }} required>
          <InputLabel>Type</InputLabel>
          <Select
            value={selectedType || ""}
            label="Type"
            onChange={handleTypeChange}
            onBlur={handleFieldBlur}
            disabled={isLocked}
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
