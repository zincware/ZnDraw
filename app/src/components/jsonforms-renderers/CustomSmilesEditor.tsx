// src/components/jsonforms-renderers/CustomSmilesEditor.tsx
// Based on https://github.com/epam/ketcher/blob/master/example/src/App.tsx

import { withJsonFormsControlProps } from "@jsonforms/react";
import { Box, Button, TextField } from "@mui/material";
import { rankWith, schemaMatches, type ControlProps } from "@jsonforms/core";
import { useState, useEffect } from "react";
import SmilesEditDialog from "./SmilesEditDialog";
import MoleculePreview from "../shared/MoleculePreview";
import FormLabelWithHelp from "../shared/FormLabelWithHelp";
import { LAYOUT_CONSTANTS } from "../../constants/layout";

/**
 * Custom SMILES editor component using Ketcher for molecular structure editing.
 */
const CustomSmilesEditor = ({
  data,
  handleChange,
  path,
  label,
  schema,
  required,
}: ControlProps) => {
  const [dialogOpen, setDialogOpen] = useState(false);
  const [currentSmiles, setCurrentSmiles] = useState(data ?? "");

  // Sync local state with incoming data prop (e.g., after form reset)
  useEffect(() => {
    setCurrentSmiles(data ?? "");
  }, [data]);

  const handleOpenDialog = () => {
    setDialogOpen(true);
  };

  const handleSaveEdit = (newSmiles: string) => {
    setCurrentSmiles(newSmiles);
    handleChange(path, newSmiles);
    setDialogOpen(false);
  };

  const handleTextChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newSmiles = event.target.value;
    setCurrentSmiles(newSmiles);
    handleChange(path, newSmiles);
  };


  return (
    <Box sx={{ marginBottom: 2 }}>
      {/* Form Label with Help Icon */}
      <FormLabelWithHelp
        label={label}
        required={required}
        helpText={schema.description}
      />

      {/* SMILES Input and Draw Button */}
      <Box sx={{ display: "flex", gap: 1, marginTop: 1 }}>
        <TextField
          value={currentSmiles}
          onChange={handleTextChange}
          placeholder="Enter SMILES notation"
          fullWidth
          size="small"
          variant="outlined"
        />
        <Button variant="outlined" onClick={handleOpenDialog}>
          Draw
        </Button>
      </Box>

      {/* Molecule Image Preview */}
      {currentSmiles && (
        <Box
          sx={{
            marginTop: 2,
            display: "flex",
            justifyContent: "center",
            alignItems: "center",
            minHeight: 100,
            "& img": {
              border: "1px solid #e0e0e0",
              borderRadius: 1,
              padding: 1,
              maxWidth: `calc(100vw - ${LAYOUT_CONSTANTS.PRIMARY_DRAWER_WIDTH}px - 48px)`,
            },
          }}
        >
          <MoleculePreview
            smiles={currentSmiles}
            size="large"
            alt="Molecule structure"
            throttleMs={500}
          />
        </Box>
      )}

      {/* Ketcher Editor Dialog */}
      <SmilesEditDialog
        open={dialogOpen}
        initialSmiles={currentSmiles}
        title="Molecular Structure Editor"
        onSave={handleSaveEdit}
        onCancel={() => setDialogOpen(false)}
      />
    </Box>
  );
};

/**
 * Tester function to determine when to use this custom SMILES editor.
 * Matches schemas with x-custom-type === "smiles" and type === "string".
 */
export const customSmilesEditorTester = rankWith(
  5, // High rank ensures it's chosen over default string renderers
  schemaMatches((schema) => {
    return (
      (schema as any)?.["x-custom-type"] === "smiles" &&
      schema.type === "string"
    );
  })
);

export default withJsonFormsControlProps(CustomSmilesEditor);
