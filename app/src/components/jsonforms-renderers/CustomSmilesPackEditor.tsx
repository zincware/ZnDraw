// src/components/jsonforms-renderers/CustomSmilesPackEditor.tsx
// Table-based editor for packing multiple molecules with SMILES notation

import { withJsonFormsControlProps } from "@jsonforms/react";
import { rankWith, schemaMatches, type ControlProps } from "@jsonforms/core";
import { useState, useEffect, useCallback, useMemo } from "react";
import {
  Box,
  Button,
  IconButton,
  Tooltip,
  Alert,
  Typography,
  Paper,
} from "@mui/material";
import AddIcon from "@mui/icons-material/Add";
import DeleteIcon from "@mui/icons-material/Delete";
import EditIcon from "@mui/icons-material/Edit";
import {
  DataGrid,
  GridColDef,
  GridRowsProp,
  GridRowModel,
} from "@mui/x-data-grid";
import SmilesEditDialog from "./SmilesEditDialog";
import MoleculePreview from "../shared/MoleculePreview";
import FormLabelWithHelp from "../shared/FormLabelWithHelp";

interface MoleculeRow {
  id: number;
  smiles: string;
  count: number;
}

/**
 * Custom SMILES pack editor for building molecular boxes.
 * Allows editing multiple molecules with their counts in a table format.
 */
const CustomSmilesPackEditor = ({
  data,
  handleChange,
  path,
  label,
  schema,
  required,
}: ControlProps) => {
  // Normalize data to array of MoleculeRow
  const [molecules, setMolecules] = useState<MoleculeRow[]>(() => {
    if (!data || !Array.isArray(data)) return [];
    return data.map((item: any, idx: number) => ({
      id: idx,
      smiles: item?.smiles || "",
      count: item?.count || 1,
    }));
  });

  const [editDialogOpen, setEditDialogOpen] = useState(false);
  const [editingMolecule, setEditingMolecule] = useState<MoleculeRow | null>(null);

  // Sync local state with incoming data prop
  useEffect(() => {
    if (!data || !Array.isArray(data)) {
      setMolecules([]);
      return;
    }
    const normalized = data.map((item: any, idx: number) => ({
      id: idx,
      smiles: item?.smiles || "",
      count: item?.count || 1,
    }));
    setMolecules(normalized);
  }, [data]);

  // Update parent when molecules change
  const updateParent = useCallback(
    (newMolecules: MoleculeRow[]) => {
      const output = newMolecules.map((m) => ({
        smiles: m.smiles,
        count: m.count,
      }));
      handleChange(path, output);
    },
    [handleChange, path]
  );

  // Add new molecule
  const handleAddMolecule = useCallback(() => {
    const newId = molecules.length > 0 ? Math.max(...molecules.map((m) => m.id)) + 1 : 0;
    const newMolecules = [
      ...molecules,
      { id: newId, smiles: "", count: 1 },
    ];
    setMolecules(newMolecules);
    updateParent(newMolecules);
  }, [molecules, updateParent]);


  // Open edit dialog
  const handleOpenEdit = useCallback((molecule: MoleculeRow) => {
    setEditingMolecule(molecule);
    setEditDialogOpen(true);
  }, []);

  // Save edited SMILES
  const handleSaveEdit = useCallback(
    (newSmiles: string) => {
      if (!editingMolecule) return;
      const newMolecules = molecules.map((m) =>
        m.id === editingMolecule.id ? { ...m, smiles: newSmiles } : m
      );
      setMolecules(newMolecules);
      updateParent(newMolecules);
      setEditDialogOpen(false);
      setEditingMolecule(null);
    },
    [editingMolecule, molecules, updateParent]
  );


  // Handle cell edit
  const processRowUpdate = useCallback(
    (newRow: GridRowModel<MoleculeRow>) => {
      const newMolecules = molecules.map((m) =>
        m.id === newRow.id
          ? {
              ...m,
              smiles: newRow.smiles,
              count: Number(newRow.count) || 1,
            }
          : m
      );
      setMolecules(newMolecules);
      updateParent(newMolecules);
              return newRow;
    },
    [molecules, updateParent]
  );

  // Validation
  const validationErrors = useMemo(() => {
    const errors: string[] = [];
    molecules.forEach((m, idx) => {
      if (!m.smiles || m.smiles.trim() === "") {
        errors.push(`Row ${idx + 1}: SMILES is required`);
      }
      if (m.count < 1) {
        errors.push(`Row ${idx + 1}: Count must be at least 1`);
      }
    });
    return errors;
  }, [molecules]);

  // Define columns
  const columns: GridColDef<MoleculeRow>[] = useMemo(
    () => [
      {
        field: "count",
        headerName: "Count",
        width: 80,
        editable: true,
        type: "number",
      },
      {
        field: "smiles",
        headerName: "SMILES",
        flex: 1,
        minWidth: 150,
        editable: true,
        type: "string",
        renderCell: (params) => {
          if (!params.value || params.value.trim() === "") {
            return (
              <Typography
                variant="body2"
                sx={{
                  color: "text.disabled",
                  fontStyle: "italic",
                }}
              >
                Click to enter SMILES
              </Typography>
            );
          }
          return params.value;
        },
      },
      {
        field: "preview",
        headerName: "Preview",
        width: 120,
        sortable: false,
        filterable: false,
        editable: false,
        renderCell: (params) => (
          <MoleculePreview
            smiles={params.row.smiles}
            size="small"
            onClick={() => handleOpenEdit(params.row)}
            errorAsTooltip={true}
          />
        ),
      },
      {
        field: "actions",
        headerName: "Actions",
        width: 120,
        sortable: false,
        filterable: false,
        editable: false,
        renderCell: (params) => (
          <Box sx={{ display: "flex", gap: 0.5 }}>
            <Tooltip title="Open molecular editor">
              <IconButton
                size="small"
                onClick={() => handleOpenEdit(params.row)}
                color="primary"
              >
                <EditIcon fontSize="small" />
              </IconButton>
            </Tooltip>
            <Tooltip title="Delete molecule">
              <IconButton
                size="small"
                onClick={() => {
                  const newMolecules = molecules.filter((m) => m.id !== params.row.id);
                  setMolecules(newMolecules);
                  updateParent(newMolecules);
                }}
                color="error"
              >
                <DeleteIcon fontSize="small" />
              </IconButton>
            </Tooltip>
          </Box>
        ),
      },
    ],
    [molecules, handleOpenEdit, updateParent]
  );

  const rows: GridRowsProp<MoleculeRow> = molecules;

  return (
    <Box sx={{ marginBottom: 2 }}>
      {/* Header */}
      <FormLabelWithHelp
        label={label}
        required={required}
        helpText={schema.description}
      />

      {/* DataGrid Table */}
      <Paper sx={{ mb: 2 }}>
        <DataGrid
          rows={rows}
          columns={columns}
          processRowUpdate={processRowUpdate}
          onProcessRowUpdateError={(error) => {
            console.error("Row update error:", error);
          }}
          getRowHeight={() => 80}
          getCellClassName={(params) => {
            // Highlight invalid SMILES cells
            if (params.field === "smiles") {
              if (!params.value || params.value.trim() === "") {
                return "invalid-cell";
              }
            }
            // Highlight invalid count cells
            if (params.field === "count") {
              if (params.value < 1) {
                return "invalid-cell";
              }
            }
            return "";
          }}
          hideFooter
          slots={{
            noRowsOverlay: () => (
            <Box
                sx={{
                  display: "flex",
                  flexDirection: "column",
                  alignItems: "center",
                  justifyContent: "center",
                  height: "100%",
                  gap: 2,
                }}
              >
                <Typography variant="body2" color="text.secondary">
                  No molecules added yet
                </Typography>
                <Button
                  size="small"
                  variant="outlined"
                  startIcon={<AddIcon />}
                  onClick={handleAddMolecule}
                >
                  Add Molecule
                </Button>
              </Box>
            ),
          }}
          sx={{
            "& .MuiDataGrid-cell:focus": {
              outline: "none",
            },
            "& .MuiDataGrid-cell:focus-within": {
              outline: "2px solid #1976d2",
            },
            "& .MuiDataGrid-cell": {
              display: "flex",
              alignItems: "center",
            },
            "& .invalid-cell": {
              border: "1px solid #d32f2f",
              backgroundColor: "rgba(211, 47, 47, 0.05)",
            },
          }}
        />
        {/* Add Molecule Button inside Paper - only show if molecules exist */}
        {molecules.length > 0 && (
          <Box
            sx={{
              borderTop: "1px solid",
              borderColor: "divider",
              p: 1.5,
              display: "flex",
              justifyContent: "center",
            }}
          >
            <Button
              size="small"
              startIcon={<AddIcon />}
              onClick={handleAddMolecule}
              variant="outlined"
            >
              Add Molecule
            </Button>
          </Box>
        )}
      </Paper>


      {/* Edit Dialog */}
      <SmilesEditDialog
        open={editDialogOpen}
        initialSmiles={editingMolecule?.smiles || ""}
        title={`Edit Molecule ${editingMolecule ? molecules.findIndex((m) => m.id === editingMolecule.id) + 1 : ""}`}
        onSave={handleSaveEdit}
        onCancel={() => {
          setEditDialogOpen(false);
          setEditingMolecule(null);
        }}
      />
    </Box>
  );
};

/**
 * Tester function to determine when to use this custom SMILES pack editor.
 * Matches schemas with x-custom-type === "smiles-pack" or "smiles-pack-array".
 */
export const customSmilesPackEditorTester = rankWith(
  6, // Higher priority than default array renderer
  schemaMatches((schema) => {
    return (
      (schema as any)?.["x-custom-type"] === "smiles-pack-array" &&
      schema.type === "array"
    );
  })
);

export default withJsonFormsControlProps(CustomSmilesPackEditor);
