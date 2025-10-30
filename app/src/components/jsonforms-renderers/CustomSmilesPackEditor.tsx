// src/components/jsonforms-renderers/CustomSmilesPackEditor.tsx
// Table-based editor for packing multiple molecules with SMILES notation

import { withJsonFormsControlProps } from "@jsonforms/react";
import { rankWith, schemaMatches, type ControlProps } from "@jsonforms/core";
import { useState, useEffect, useCallback, useMemo } from "react";
import {
  Box,
  Button,
  FormLabel,
  IconButton,
  Tooltip,
  Alert,
  Typography,
  Chip,
  CircularProgress,
  Paper,
} from "@mui/material";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";
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
import { convertMoleculeToImage } from "../../myapi/client";

interface MoleculeRow {
  id: number;
  smiles: string;
  count: number;
}

interface MoleculeImage {
  smiles: string;
  image: string | null;
  loading: boolean;
  error: string | null;
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
  const [moleculeImages, setMoleculeImages] = useState<Map<string, MoleculeImage>>(
    new Map()
  );

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

  // Fetch molecule image
  const fetchMoleculeImage = useCallback(async (smiles: string) => {
    if (!smiles || smiles.trim() === "") return;

    setMoleculeImages((prev) => {
      const updated = new Map(prev);
      updated.set(smiles, {
        smiles,
        image: null,
        loading: true,
        error: null,
      });
      return updated;
    });

    try {
      const response = await convertMoleculeToImage({
        type: "smiles",
        data: smiles,
      });
      setMoleculeImages((prev) => {
        const updated = new Map(prev);
        updated.set(smiles, {
          smiles,
          image: response.image,
          loading: false,
          error: null,
        });
        return updated;
      });
    } catch (error: any) {
      console.error("Error fetching molecule image:", error);
      setMoleculeImages((prev) => {
        const updated = new Map(prev);
        updated.set(smiles, {
          smiles,
          image: null,
          loading: false,
          error: error.response?.data?.error || "Failed to generate image",
        });
        return updated;
      });
    }
  }, []);

  // Load images automatically for all molecules with SMILES
  useEffect(() => {
    molecules.forEach((m) => {
      if (m.smiles && m.smiles.trim() !== "" && !moleculeImages.has(m.smiles)) {
        fetchMoleculeImage(m.smiles);
      }
    });
  }, [molecules, moleculeImages, fetchMoleculeImage]);

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

  // Calculate totals
  const totalMolecules = molecules.reduce((sum, m) => sum + m.count, 0);

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
      },
      {
        field: "preview",
        headerName: "Preview",
        width: 120,
        sortable: false,
        filterable: false,
        editable: false,
        renderCell: (params) => {
          const imageData = moleculeImages.get(params.row.smiles);
          if (!params.row.smiles || params.row.smiles.trim() === "") {
            return null;
          }
          if (imageData?.loading) {
            return <CircularProgress size={24} />;
          }
          if (imageData?.error) {
            return (
              <Tooltip title={imageData.error}>
                <Typography variant="caption" color="error">
                  Error
                </Typography>
              </Tooltip>
            );
          }
          if (imageData?.image) {
            return (
              <Box
                component="img"
                src={imageData.image}
                alt={params.row.smiles}
                onClick={() => handleOpenEdit(params.row)}
                sx={{
                  width: "100%",
                  height: "60px",
                  objectFit: "contain",
                  backgroundColor: "white",
                  cursor: "pointer",
                  "&:hover": {
                    opacity: 0.8,
                  },
                }}
              />
            );
          }
          return null;
        },
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
    [molecules, handleOpenEdit, moleculeImages]
  );

  const rows: GridRowsProp<MoleculeRow> = molecules;

  return (
    <Box sx={{ marginBottom: 2 }}>
      {/* Header */}
      <Box sx={{ display: "flex", alignItems: "center", gap: 1, mb: 1 }}>
        <FormLabel required={required}>{label}</FormLabel>
        {schema.description && (
          <Tooltip title={schema.description} placement="top">
            <IconButton size="small" sx={{ padding: 0.5 }}>
              <HelpOutlineIcon fontSize="small" />
            </IconButton>
          </Tooltip>
        )}
      </Box>

      {/* Validation Errors */}
      {validationErrors.length > 0 && (
        <Alert severity="error" sx={{ mb: 2 }}>
          {validationErrors.map((err, idx) => (
            <div key={idx}>{err}</div>
          ))}
        </Alert>
      )}

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
