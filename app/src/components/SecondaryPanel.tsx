import { useEffect, useMemo, useState, useRef, useCallback, memo } from "react";
import {
  Box,
  Typography,
  Divider,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Button,
  SelectChangeEvent,
  Fade,
  Chip,
} from "@mui/material";
import SaveIcon from "@mui/icons-material/Save";
import { JsonForms } from "@jsonforms/react";
import { materialCells } from "@jsonforms/material-renderers";
import { useFormStore } from "../formStore";
import {
  useSchemas,
  useExtensionData,
  useSubmitExtension,
  useFrameKeys,
  useExtensionStats,
} from "../hooks/useSchemas";
import { useAppStore } from "../store";
import { ExtensionStatusChips } from "./ExtensionStatusChips";
import { JobHistoryPanel } from "./JobHistoryPanel";
import { debounce } from "lodash";
import { customRenderers, injectDynamicEnums, schemaRequiresMetadata } from "../utils/jsonforms";
import { PanelSkeleton, FormSkeleton } from "./shared/LoadingSkeletons";

interface SecondaryPanelProps {
  panelTitle: string;
}

const SecondaryPanel = ({ panelTitle }: SecondaryPanelProps) => {
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const userName = useAppStore((state) => state.userName);
  const geometries = useAppStore((state) => state.geometries);
  const [localFormData, setLocalFormData] = useState<any>({});
  const ignoreFirstChangeRef = useRef(true);

  if (!roomId || !userName) {
    return <Typography sx={{ p: 2 }}>Joining room...</Typography>;
  }

  const { selectedExtensions, setSelectedExtension } = useFormStore();
  const selectedExtensionKey = selectedExtensions[panelTitle] || null;

  const {
    data: schemas,
    isLoading: isLoadingSchemas,
    isError: isSchemasError,
  } = useSchemas(roomId, panelTitle);

  // Parse composite key "name:public" to extract name and public flag
  const parseExtensionKey = (key: string | null): { name: string; public: boolean } | null => {
    if (!key) return null;
    const [name, publicStr] = key.split(":");
    return { name, public: publicStr === "true" };
  };

  const selectedExtension = parseExtensionKey(selectedExtensionKey);
  const selectedExtensionName = selectedExtension?.name || null;

  // Find the selected extension object in the schemas list
  const selectedExtensionObject = useMemo(() => {
    if (!schemas || !selectedExtension) return null;
    return schemas.find(
      ext => ext.name === selectedExtension.name && ext.public === selectedExtension.public
    );
  }, [schemas, selectedExtension]);

  // Check if the selected extension's schema requires metadata
  const currentSchema = selectedExtensionObject?.schema;
  const needsMetadata = useMemo(
    () => currentSchema ? schemaRequiresMetadata(currentSchema) : false,
    [currentSchema]
  );

  // --- MODIFICATION: Fetch the frame keys only when needed (deferred loading) ---
  const { data: metadata, isLoading: isLoadingMetadata } = useFrameKeys(
    roomId,
    0,
    needsMetadata // Only fetch if the current schema needs it
  );

  const {
    data: serverData,
    isLoading: isLoadingData,
    isError: isDataError,
  } = useExtensionData(roomId, panelTitle, selectedExtensionName || "");

  // Fetch extension stats for worker availability and job counts
  const scope = selectedExtension?.public ? "global" : "user";
  const { data: extensionStats } = useExtensionStats(
    roomId,
    scope,
    panelTitle,
    selectedExtensionName || ""
  );

  useEffect(() => {
    if (!isLoadingData && serverData !== undefined) {
      setLocalFormData(serverData ?? {});
      ignoreFirstChangeRef.current = true;
    }
  }, [isLoadingData, serverData, selectedExtensionKey]);

  const { mutate: submit, isPending: isSubmitting } = useSubmitExtension();

  const debouncedSubmit = useMemo(
    () =>
      debounce((data: any) => {
        if (!selectedExtensionName || !selectedExtension || !roomId || !userName || !schemas) return;
        submit({
          roomId,
          userName,
          category: panelTitle,
          extension: selectedExtensionName,
          data: data,
          isPublic: selectedExtension.public,
        });
      }, 100), // Increased debounce time slightly
    [selectedExtensionName, selectedExtension, roomId, userName, panelTitle, submit, schemas],
  );

  useEffect(() => {
    return () => {
      debouncedSubmit.cancel();
    };
  }, [debouncedSubmit]);

  const handleFormChange = useCallback(
    ({ data }: { data: any }) => {
      const safeData = data ?? {};
      if (ignoreFirstChangeRef.current) {
        ignoreFirstChangeRef.current = false;
        return;
      }
      setLocalFormData(safeData);

      if (panelTitle === "settings") {
        debouncedSubmit(safeData);
      }
    },
    [panelTitle, debouncedSubmit],
  );

  const handleSelectionChange = (event: SelectChangeEvent<string>) => {
    setSelectedExtension(panelTitle, event.target.value || null);
  };

  const handleSubmit = () => {
    if (!selectedExtensionName || !selectedExtension || !roomId || !userName) return;

    submit({
      roomId,
      category: panelTitle,
      extension: selectedExtensionName,
      data: localFormData,
      isPublic: selectedExtension.public,
    });
  };

  // --- MODIFICATION: Create a dynamic schema by injecting metadata ---
  const dynamicSchema = useMemo(() => {
    const originalSchema = selectedExtensionObject?.schema;
    if (!originalSchema) return null;

    // Call our new helper function to handle injection generically
    return injectDynamicEnums(originalSchema, metadata, geometries);

  }, [selectedExtensionObject, metadata, geometries]);

  // Create form options with composite keys and check for duplicates
  const formOptions = useMemo(() => {
    if (!schemas || !Array.isArray(schemas)) return [];

    // Check if any name appears with both public=true and public=false
    const nameCounts = new Map<string, { hasPublic: boolean; hasPrivate: boolean }>();

    schemas.forEach(ext => {
      const existing = nameCounts.get(ext.name) || { hasPublic: false, hasPrivate: false };
      if (ext.public) {
        existing.hasPublic = true;
      } else {
        existing.hasPrivate = true;
      }
      nameCounts.set(ext.name, existing);
    });

    // Build options with composite keys
    return schemas.map(ext => {
      const counts = nameCounts.get(ext.name)!;
      const showBadge = counts.hasPublic && counts.hasPrivate;

      return {
        name: ext.name,
        compositeKey: `${ext.name}:${ext.public}`,
        showBadge,
        isPublic: ext.public,
      };
    });
  }, [schemas]);

  if (isSchemasError || isDataError) {
    return (
      <Typography color="error" sx={{ p: 2 }}>
        Failed to load configuration.
      </Typography>
    );
  }

  return (
    <Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
      <Typography variant="h6" sx={{ p: 2, pb: 1, flexShrink: 0 }}>
        {panelTitle}
      </Typography>
      <Divider />

      <Box sx={{ p: 2, pb: 12, flexGrow: 1, overflowY: "auto" }}>
        {isLoadingSchemas ? (
          <PanelSkeleton />
        ) : (
          <>
            <FormControl fullWidth sx={{ mb: 2 }}>
              <InputLabel id="panel-select-label">{panelTitle} Method</InputLabel>
              <Select
                labelId="panel-select-label"
                value={selectedExtensionKey || ""}
                label={`${panelTitle} Method`}
                onChange={handleSelectionChange}
              >
                {formOptions.map((option) => (
                  <MenuItem key={option.compositeKey} value={option.compositeKey}>
                    <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
                      <Typography>{option.name}</Typography>
                      {option.showBadge && (
                        <Chip
                          label={option.isPublic ? "Public" : "Private"}
                          size="small"
                          color={option.isPublic ? "success" : "default"}
                        />
                      )}
                    </Box>
                  </MenuItem>
                ))}
              </Select>
            </FormControl>

            {selectedExtensionObject && (
              <ExtensionStatusChips
                metadata={selectedExtensionObject}
                stats={extensionStats}
              />
            )}

            {selectedExtension && (
              <>
                {panelTitle !== "settings" && (
                  <Button
                    variant="contained"
                    startIcon={<SaveIcon />}
                    onClick={handleSubmit}
                    disabled={isSubmitting || isLoadingData || isLoadingMetadata}
                    fullWidth
                    color="primary"
                    sx={{ mb: 2 }}
                  >
                    {isSubmitting ? "Running..." : "Run Extension"}
                  </Button>
                )}

                {isLoadingData || isLoadingMetadata ? (
                  <FormSkeleton />
                ) : dynamicSchema ? (
                  <Fade in={!isLoadingData && !isLoadingMetadata} timeout={200}>
                    <Box>
                      <JsonForms
                        key={selectedExtensionKey}
                        schema={dynamicSchema}
                        data={localFormData}
                        renderers={customRenderers}
                        cells={materialCells}
                        onChange={handleFormChange}
                      />
                    </Box>
                  </Fade>
                ) : null}

                {/* Job History Panel - only show for non-settings panels */}
                {panelTitle !== "settings" && selectedExtension && (
                  <JobHistoryPanel
                    roomId={roomId}
                    category={panelTitle}
                    extension={selectedExtensionName}
                    isPublic={selectedExtension.public}
                  />
                )}
              </>
            )}
          </>
        )}
      </Box>
    </Box>
  );
};

// Memoize to prevent unnecessary re-renders
export default memo(SecondaryPanel);