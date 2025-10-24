import { useState, useMemo, useCallback } from "react";
import {
  Box,
  Chip,
  IconButton,
  TextField,
  Button,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogContentText,
  DialogActions,
  Tooltip,
  Typography,
} from "@mui/material";
import { DataGrid, GridColDef, GridRenderCellParams, GridRowParams } from "@mui/x-data-grid";
import EditIcon from "@mui/icons-material/Edit";
import DeleteIcon from "@mui/icons-material/Delete";
import ClearIcon from "@mui/icons-material/Clear";
import RadioButtonCheckedIcon from "@mui/icons-material/RadioButtonChecked";
import RadioButtonUncheckedIcon from "@mui/icons-material/RadioButtonUnchecked";
import { useAppStore } from "../store";
import {
  createUpdateSelectionGroup,
  deleteSelectionGroup,
  loadSelectionGroup,
} from "../myapi/client";

interface SelectionGroupRow {
  id: string;
  name: string;
  isCurrent: boolean;
  isActive: boolean;
  selections: Record<string, number>; // geometry -> count
  rawData: Record<string, number[]>; // geometry -> indices array
}

export default function SelectionGroupsPanel() {
  const { roomId, selections, selectionGroups, activeSelectionGroup, updateSelectionForGeometry, showSnackbar } =
    useAppStore();

  // Local state
  const [currentGroupName, setCurrentGroupName] = useState("");
  const [editingGroupName, setEditingGroupName] = useState<string | null>(null);
  const [editNameValue, setEditNameValue] = useState("");
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [groupToDelete, setGroupToDelete] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(false);

  // Get all unique geometry keys across all groups and current selection
  const allGeometryKeys = useMemo(() => {
    const keys = new Set<string>();

    // Add keys from current selection
    Object.keys(selections).forEach((key) => {
      if (selections[key] && selections[key].length > 0) {
        keys.add(key);
      }
    });

    // Add keys from all saved groups
    Object.values(selectionGroups).forEach((group) => {
      Object.keys(group).forEach((key) => {
        if (group[key] && group[key].length > 0) {
          keys.add(key);
        }
      });
    });

    return Array.from(keys).sort();
  }, [selections, selectionGroups]);

  // Prepare rows for DataGrid
  const rows = useMemo((): SelectionGroupRow[] => {
    const result: SelectionGroupRow[] = [];

    // Current selection row (always first)
    const currentSelectionCounts: Record<string, number> = {};
    allGeometryKeys.forEach((key) => {
      currentSelectionCounts[key] = selections[key]?.length || 0;
    });

    result.push({
      id: "current",
      name: "", // Will show input field
      isCurrent: true,
      isActive: false,
      selections: currentSelectionCounts,
      rawData: selections,
    });

    // Saved group rows
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

  // Generate dynamic columns
  const columns = useMemo((): GridColDef[] => {
    const cols: GridColDef[] = [
      // Name column with radio button
      {
        field: "name",
        headerName: "Name",
        width: 200,
        flex: 1,
        renderCell: (params: GridRenderCellParams<SelectionGroupRow>) => {
          const row = params.row;

          if (row.isCurrent) {
            // Current selection row - show input
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

          // Saved group row
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

    // Add dynamic geometry columns
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

    // Actions column
    cols.push({
      field: "actions",
      headerName: "Actions",
      width: 100,
      sortable: false,
      renderCell: (params: GridRenderCellParams<SelectionGroupRow>) => {
        const row = params.row;

        if (row.isCurrent) {
          // Current selection - Clear button
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

        // Saved group - Edit and Delete buttons
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

  // Check if current selection has any data
  const hasCurrentSelection = useMemo(
    () => Object.values(selections).some((indices) => indices && indices.length > 0),
    [selections]
  );

  // Handlers
  const handleSaveCurrentSelection = useCallback(async () => {
    if (!roomId || !currentGroupName.trim()) return;

    // Check for duplicate name
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

      // Check for duplicate name
      if (selectionGroups[editNameValue.trim()] && editNameValue.trim() !== oldName) {
        setSnackbar({
          open: true,
          message: `Group "${editNameValue.trim()}" already exists`,
          severity: "error",
        });
        setEditingGroupName(null);
        return;
      }

      setIsLoading(true);
      try {
        // Delete old group
        await deleteSelectionGroup(roomId, oldName);
        // Create new group with new name
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

      // Don't load current selection row or if already editing
      if (row.isCurrent || editingGroupName) return;

      if (!roomId) return;

      setIsLoading(true);
      try {
        await loadSelectionGroup(roomId, row.id);
        showSnackbar(`Loaded group "${row.name}"`, "success");
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

    // Clear all geometry selections
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

  return (
    <Box sx={{ height: "100%", width: "100%", display: "flex", flexDirection: "column", p: 2 }}>
      <Typography variant="h6" sx={{ mb: 1 }}>
        Selection Groups
      </Typography>

      {/* DataGrid */}
      <Box sx={{ flexGrow: 1, minHeight: 0 }}>
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

    </Box>
  );
}
