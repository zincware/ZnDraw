import { useEffect, useMemo, useState, useRef, useCallback } from "react";
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
  Chip,
  IconButton,
  TextField,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogContentText,
  DialogActions,
  Tooltip,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Fade,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import SaveIcon from "@mui/icons-material/Save";
import EditIcon from "@mui/icons-material/Edit";
import DeleteIcon from "@mui/icons-material/Delete";
import ClearIcon from "@mui/icons-material/Clear";
import RadioButtonCheckedIcon from "@mui/icons-material/RadioButtonChecked";
import RadioButtonUncheckedIcon from "@mui/icons-material/RadioButtonUnchecked";
import { DataGrid, GridColDef, GridRenderCellParams, GridRowParams } from "@mui/x-data-grid";
import { JsonForms } from "@jsonforms/react";
import { materialCells } from "@jsonforms/material-renderers";
import { debounce } from "lodash";
import { useFormStore } from "../formStore";
import { useAppStore } from "../store";
import {
  useSchemas,
  useExtensionData,
  useSubmitExtension,
  useFrameMetadata,
} from "../hooks/useSchemas";
import {
  createUpdateSelectionGroup,
  deleteSelectionGroup,
  loadSelectionGroup,
} from "../myapi/client";
import { ExtensionStatusChips } from "./ExtensionStatusChips";
import { customRenderers, injectDynamicEnums, schemaRequiresMetadata } from "../utils/jsonforms";
import { SelectionToolsSkeleton, FormSkeleton } from "./shared/LoadingSkeletons";

interface SelectionGroupRow {
  id: string;
  name: string;
  isCurrent: boolean;
  isActive: boolean;
  selections: Record<string, number>; // geometry -> count
  rawData: Record<string, number[]>; // geometry -> indices array
}

export default function SelectionsPanel() {
  const { roomId, userName, selections, selectionGroups, activeSelectionGroup, updateSelectionForGeometry, geometries, showSnackbar } =
    useAppStore();

  // Selection Groups state
  const [currentGroupName, setCurrentGroupName] = useState("");
  const [editingGroupName, setEditingGroupName] = useState<string | null>(null);
  const [editNameValue, setEditNameValue] = useState("");
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [groupToDelete, setGroupToDelete] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(false);

  // Selection Tools state
  const [localFormData, setLocalFormData] = useState<any>({});
  const ignoreFirstChangeRef = useRef(true);
  const scrollContainerRef = useRef<HTMLDivElement>(null);
  const panelTitle = "selections";

  const { selectedExtensions, setSelectedExtension } = useFormStore();
  const selectedExtension = selectedExtensions[panelTitle] || null;

  const {
    data: schemas,
    isLoading: isLoadingSchemas,
    isError: isSchemasError,
  } = useSchemas(roomId!, panelTitle);

  // Check if the selected extension's schema requires metadata
  const currentSchema = schemas?.[selectedExtension ?? ""]?.schema;
  const needsMetadata = useMemo(
    () => currentSchema ? schemaRequiresMetadata(currentSchema) : false,
    [currentSchema]
  );

  const { data: metadata, isLoading: isLoadingMetadata } = useFrameMetadata(
    roomId!,
    0,
    needsMetadata
  );

  const {
    data: serverData,
    isLoading: isLoadingData,
    isError: isDataError,
  } = useExtensionData(roomId!, userName!, panelTitle, selectedExtension || "");

  useEffect(() => {
    if (!isLoadingData && serverData !== undefined) {
      setLocalFormData(serverData ?? {});
      ignoreFirstChangeRef.current = true;
    }
  }, [isLoadingData, serverData, selectedExtension]);

  const { mutate: submit, isPending: isSubmitting } = useSubmitExtension();

  const debouncedSubmit = useMemo(
    () =>
      debounce((data: any) => {
        if (!selectedExtension || !roomId || !userName) return;
        submit({
          roomId,
          userName,
          category: panelTitle,
          extension: selectedExtension,
          data: data,
        });
      }, 100),
    [selectedExtension, roomId, userName, submit],
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
    },
    [debouncedSubmit],
  );

  const handleSelectionChange = (event: SelectChangeEvent<string>) => {
    setSelectedExtension(panelTitle, event.target.value || null);
  };

  const handleSubmit = () => {
    if (!selectedExtension || !roomId || !userName) return;
    submit({
      roomId,
      userName,
      category: panelTitle,
      extension: selectedExtension,
      data: localFormData,
    });
  };

  const dynamicSchema = useMemo(() => {
    const originalSchema = schemas?.[selectedExtension ?? ""]?.schema;
    if (!originalSchema) return null;
    return injectDynamicEnums(originalSchema, metadata, geometries);
  }, [schemas, selectedExtension, metadata, geometries]);

  const formOptions = useMemo(() => Object.keys(schemas || {}), [schemas]);

  // Selection Groups logic
  const allGeometryKeys = useMemo(() => {
    const keys = new Set<string>();

    Object.keys(selections).forEach((key) => {
      if (selections[key] && selections[key].length > 0) {
        keys.add(key);
      }
    });

    Object.values(selectionGroups).forEach((group) => {
      Object.keys(group).forEach((key) => {
        if (group[key] && group[key].length > 0) {
          keys.add(key);
        }
      });
    });

    return Array.from(keys).sort();
  }, [selections, selectionGroups]);

  const rows = useMemo((): SelectionGroupRow[] => {
    const result: SelectionGroupRow[] = [];

    const currentSelectionCounts: Record<string, number> = {};
    allGeometryKeys.forEach((key) => {
      currentSelectionCounts[key] = selections[key]?.length || 0;
    });

    result.push({
      id: "current",
      name: "",
      isCurrent: true,
      isActive: false,
      selections: currentSelectionCounts,
      rawData: selections,
    });

    const groupNames = Object.keys(selectionGroups).sort();
    groupNames.forEach((groupName) => {
      const groupData = selectionGroups[groupName];
      const groupCounts: Record<string, number> = {};

      allGeometryKeys.forEach((key) => {
        groupCounts[key] = groupData[key]?.length || 0;
      });

      result.push({
        id: groupName,
        name: groupName,
        isCurrent: false,
        isActive: activeSelectionGroup === groupName,
        selections: groupCounts,
        rawData: groupData,
      });
    });

    return result;
  }, [selections, selectionGroups, activeSelectionGroup, allGeometryKeys]);

  const columns = useMemo((): GridColDef[] => {
    const cols: GridColDef[] = [
      {
        field: "name",
        headerName: "Name",
        width: 200,
        flex: 1,
        renderCell: (params: GridRenderCellParams<SelectionGroupRow>) => {
          const row = params.row;

          if (row.isCurrent) {
            return (
              <Box sx={{ display: "flex", alignItems: "center", gap: 1, width: "100%", height: "100%" }}>
                <TextField
                  size="small"
                  fullWidth
                  placeholder="Type name to save..."
                  value={currentGroupName}
                  onChange={(e) => setCurrentGroupName(e.target.value)}
                  onKeyDown={(e) => {
                    if (e.key === "Enter" && currentGroupName.trim()) {
                      handleSaveCurrentSelection();
                    }
                  }}
                  onBlur={() => {
                    if (currentGroupName.trim()) {
                      handleSaveCurrentSelection();
                    }
                  }}
                  disabled={isLoading || !hasCurrentSelection}
                  sx={{ "& .MuiInputBase-root": { fontSize: "0.875rem" } }}
                />
              </Box>
            );
          }

          const isEditing = editingGroupName === row.id;

          return (
            <Box sx={{ display: "flex", alignItems: "center", gap: 1, width: "100%", height: "100%" }}>
              {row.isActive ? (
                <RadioButtonCheckedIcon sx={{ color: "success.main" }} fontSize="small" />
              ) : (
                <RadioButtonUncheckedIcon sx={{ color: "action.disabled" }} fontSize="small" />
              )}
              {isEditing ? (
                <TextField
                  size="small"
                  fullWidth
                  autoFocus
                  value={editNameValue}
                  onChange={(e) => setEditNameValue(e.target.value)}
                  onKeyDown={(e) => {
                    if (e.key === "Enter") {
                      handleSaveNameEdit(row.id);
                    } else if (e.key === "Escape") {
                      setEditingGroupName(null);
                    }
                  }}
                  onBlur={() => handleSaveNameEdit(row.id)}
                  sx={{ "& .MuiInputBase-root": { fontSize: "0.875rem" } }}
                />
              ) : (
                <Typography
                  sx={{
                    fontWeight: row.isActive ? "bold" : "normal",
                    fontSize: "0.875rem",
                    flex: 1,
                  }}
                  onDoubleClick={() => {
                    setEditingGroupName(row.id);
                    setEditNameValue(row.name);
                  }}
                >
                  {row.name}
                </Typography>
              )}
            </Box>
          );
        },
      },
    ];

    allGeometryKeys.forEach((geometryKey) => {
      cols.push({
        field: geometryKey,
        headerName: geometryKey,
        width: 100,
        align: "center",
        headerAlign: "center",
        renderCell: (params: GridRenderCellParams<SelectionGroupRow>) => {
          const row = params.row;
          const count = row.selections[geometryKey] || 0;

          if (count === 0) {
            return (
              <Box sx={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100%" }}>
                <Typography sx={{ color: "text.disabled", fontSize: "0.875rem" }}>-</Typography>
              </Box>
            );
          }

          let chipColor: "primary" | "success" | "default" = "default";
          if (row.isCurrent) {
            chipColor = "primary";
          } else if (row.isActive) {
            chipColor = "success";
          }

          return (
            <Box sx={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100%" }}>
              <Chip label={count} size="small" color={chipColor} />
            </Box>
          );
        },
      });
    });

    cols.push({
      field: "actions",
      headerName: "Actions",
      width: 100,
      sortable: false,
      renderCell: (params: GridRenderCellParams<SelectionGroupRow>) => {
        const row = params.row;

        if (row.isCurrent) {
          return (
            <Box sx={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100%" }}>
              <Tooltip title="Clear all current selections">
                <span>
                  <IconButton
                    size="small"
                    onClick={(e) => {
                      e.stopPropagation();
                      handleClearCurrent();
                    }}
                    disabled={!hasCurrentSelection || isLoading}
                  >
                    <ClearIcon fontSize="small" />
                  </IconButton>
                </span>
              </Tooltip>
            </Box>
          );
        }

        return (
          <Box sx={{ display: "flex", alignItems: "center", justifyContent: "center", gap: 0.5, height: "100%" }}>
            <Tooltip title="Edit name">
              <IconButton
                size="small"
                onClick={(e) => {
                  e.stopPropagation();
                  setEditingGroupName(row.id);
                  setEditNameValue(row.name);
                }}
                disabled={isLoading}
              >
                <EditIcon fontSize="small" />
              </IconButton>
            </Tooltip>
            <Tooltip title="Delete group">
              <IconButton
                size="small"
                onClick={(e) => {
                  e.stopPropagation();
                  setGroupToDelete(row.id);
                  setDeleteDialogOpen(true);
                }}
                disabled={isLoading}
              >
                <DeleteIcon fontSize="small" />
              </IconButton>
            </Tooltip>
          </Box>
        );
      },
    });

    return cols;
  }, [allGeometryKeys, currentGroupName, editingGroupName, editNameValue, isLoading]);

  const hasCurrentSelection = useMemo(
    () => Object.values(selections).some((indices) => indices && indices.length > 0),
    [selections]
  );

  const handleSaveCurrentSelection = useCallback(async () => {
    if (!roomId || !currentGroupName.trim()) return;

    if (selectionGroups[currentGroupName.trim()]) {
      showSnackbar(`Group "${currentGroupName.trim()}" already exists`, "error");
      return;
    }

    setIsLoading(true);
    try {
      await createUpdateSelectionGroup(roomId, currentGroupName.trim(), selections);
      setCurrentGroupName("");
      showSnackbar(`Group "${currentGroupName.trim()}" saved`, "success");
    } catch (err) {
      console.error("Failed to save selection group:", err);
      showSnackbar("Failed to save selection group", "error");
    } finally {
      setIsLoading(false);
    }
  }, [roomId, currentGroupName, selections, selectionGroups, showSnackbar]);

  const handleSaveNameEdit = useCallback(
    async (oldName: string) => {
      if (!roomId || !editNameValue.trim() || editNameValue.trim() === oldName) {
        setEditingGroupName(null);
        return;
      }

      if (selectionGroups[editNameValue.trim()] && editNameValue.trim() !== oldName) {
        showSnackbar(`Group "${editNameValue.trim()}" already exists`, "error");
        setEditingGroupName(null);
        return;
      }

      setIsLoading(true);
      try {
        await deleteSelectionGroup(roomId, oldName);
        await createUpdateSelectionGroup(roomId, editNameValue.trim(), selectionGroups[oldName]);
        setEditingGroupName(null);
        showSnackbar(`Group renamed to "${editNameValue.trim()}"`, "success");
      } catch (err) {
        console.error("Failed to rename group:", err);
        showSnackbar("Failed to rename group", "error");
        setEditingGroupName(null);
      } finally {
        setIsLoading(false);
      }
    },
    [roomId, editNameValue, selectionGroups, showSnackbar]
  );

  const handleRowClick = useCallback(
    async (params: GridRowParams<SelectionGroupRow>) => {
      const row = params.row;

      if (row.isCurrent || editingGroupName) return;

      if (!roomId) return;

      setIsLoading(true);
      try {
        await loadSelectionGroup(roomId, row.id);
        // Success - no notification needed
      } catch (err) {
        console.error("Failed to load selection group:", err);
        showSnackbar(`Failed to load group "${row.name}"`, "error");
      } finally {
        setIsLoading(false);
      }
    },
    [roomId, editingGroupName, showSnackbar]
  );

  const handleClearCurrent = useCallback(() => {
    if (!hasCurrentSelection) return;

    allGeometryKeys.forEach((key) => {
      if (selections[key] && selections[key].length > 0) {
        updateSelectionForGeometry(key, []);
      }
    });
  }, [hasCurrentSelection, allGeometryKeys, selections, updateSelectionForGeometry]);

  const handleDeleteConfirm = useCallback(async () => {
    if (!roomId || !groupToDelete) return;

    setIsLoading(true);
    setDeleteDialogOpen(false);

    try {
      await deleteSelectionGroup(roomId, groupToDelete);
      showSnackbar(`Group "${groupToDelete}" deleted`, "success");
    } catch (err) {
      console.error("Failed to delete group:", err);
      showSnackbar(`Failed to delete group "${groupToDelete}"`, "error");
    } finally {
      setIsLoading(false);
      setGroupToDelete(null);
    }
  }, [roomId, groupToDelete, showSnackbar]);

  if (!roomId || !userName) {
    return <Typography sx={{ p: 2 }}>Joining room...</Typography>;
  }

  if (isSchemasError || isDataError) {
    return (
      <Typography color="error" sx={{ p: 2 }}>
        Failed to load configuration.
      </Typography>
    );
  }

  return (
    <Box sx={{ height: "100%", width: "100%", display: "flex", flexDirection: "column" }}>
      <Typography variant="h6" sx={{ p: 2, pb: 1 }}>
        Selections
      </Typography>
      <Divider />

      <Box ref={scrollContainerRef} sx={{ flexGrow: 1, overflowY: "auto", overscrollBehavior: "contain", p: 2 }}>
        {isLoadingSchemas ? (
          <SelectionToolsSkeleton />
        ) : (
          <>
        {/* Selection Tools Section */}
        <Accordion defaultExpanded>
          <AccordionSummary expandIcon={<ExpandMoreIcon />}>
            <Typography variant="subtitle1" fontWeight="medium">
              Selection Tools
            </Typography>
          </AccordionSummary>
          <AccordionDetails>
            <FormControl fullWidth sx={{ mb: 2 }}>
              <InputLabel id="panel-select-label">Selection Method</InputLabel>
              <Select
                labelId="panel-select-label"
                value={selectedExtension || ""}
                label="Selection Method"
                onChange={handleSelectionChange}
              >
                {formOptions.map((item) => (
                  <MenuItem key={item} value={item}>
                    {item}
                  </MenuItem>
                ))}
              </Select>
            </FormControl>

            {selectedExtension && schemas && schemas[selectedExtension] && (
              <ExtensionStatusChips metadata={schemas[selectedExtension]} />
            )}

            {dynamicSchema && (
              <>
                <Button
                  variant="contained"
                  startIcon={<SaveIcon />}
                  onClick={handleSubmit}
                  disabled={isSubmitting || isLoadingData}
                  fullWidth
                  color="primary"
                  sx={{ mb: 2 }}
                >
                  {isSubmitting ? "Running..." : "Run Extension"}
                </Button>

                {isLoadingData || isLoadingMetadata ? (
                  <FormSkeleton />
                ) : (
                  <Fade in={!isLoadingData && !isLoadingMetadata} timeout={200}>
                    <Box>
                      <JsonForms
                        key={selectedExtension}
                        schema={dynamicSchema}
                        data={localFormData}
                        renderers={customRenderers}
                        cells={materialCells}
                        onChange={handleFormChange}
                      />
                    </Box>
                  </Fade>
                )}
              </>
            )}
          </AccordionDetails>
        </Accordion>

        {/* Selection Groups Section */}
        <Accordion defaultExpanded sx={{ mt: 2 }}>
          <AccordionSummary expandIcon={<ExpandMoreIcon />}>
            <Typography variant="subtitle1" fontWeight="medium">
              Selection Groups
            </Typography>
          </AccordionSummary>
          <AccordionDetails sx={{ p: 0 }}>
            <Box sx={{ height: 400, width: "100%" }}>
              <DataGrid
                rows={rows}
                columns={columns}
                disableRowSelectionOnClick
                disableColumnMenu
                hideFooter
                onRowClick={handleRowClick}
                getRowClassName={(params) => {
                  if (params.row.isCurrent) return "current-selection-row";
                  if (params.row.isActive) return "active-group-row";
                  return "";
                }}
                sx={{
                  "& .MuiDataGrid-columnHeaders": {
                    minHeight: "48px !important",
                    maxHeight: "48px !important",
                  },
                  "& .MuiDataGrid-columnHeader": {
                    padding: "8px",
                  },
                  "& .current-selection-row": {
                    bgcolor: "rgba(33, 150, 243, 0.08)",
                    borderTop: "2px solid",
                    borderTopColor: "primary.main",
                    "&:hover": {
                      bgcolor: "rgba(33, 150, 243, 0.12)",
                    },
                  },
                  "& .active-group-row": {
                    bgcolor: "rgba(76, 175, 80, 0.08)",
                    "&:hover": {
                      bgcolor: "rgba(76, 175, 80, 0.12)",
                    },
                  },
                  "& .MuiDataGrid-row": {
                    cursor: "pointer",
                  },
                  "& .MuiDataGrid-row.current-selection-row": {
                    cursor: "default",
                  },
                }}
              />
            </Box>
          </AccordionDetails>
        </Accordion>

        {/* Delete Confirmation Dialog */}
        <Dialog open={deleteDialogOpen} onClose={() => setDeleteDialogOpen(false)}>
          <DialogTitle>Delete Selection Group</DialogTitle>
          <DialogContent>
            <DialogContentText>
              Are you sure you want to delete the group "{groupToDelete}"? This action cannot be
              undone.
            </DialogContentText>
          </DialogContent>
          <DialogActions>
            <Button onClick={() => setDeleteDialogOpen(false)}>Cancel</Button>
            <Button onClick={handleDeleteConfirm} color="error" variant="contained">
              Delete
            </Button>
          </DialogActions>
        </Dialog>
        </>
        )}
      </Box>
    </Box>
  );
}
