// src/components/jsonforms-renderers/CustomSmilesEditor.tsx
// Based on https://github.com/epam/ketcher/blob/master/example/src/App.tsx

import { withJsonFormsControlProps } from "@jsonforms/react";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";
import {
  Box,
  Button,
  FormLabel,
  IconButton,
  TextField,
  Tooltip,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
} from "@mui/material";
import { rankWith, schemaMatches, type ControlProps } from "@jsonforms/core";
import { StrictMode, useState, useEffect } from "react";
import { Editor, InfoModal } from "ketcher-react";
import { StructServiceProvider } from "ketcher-core";
import { StandaloneStructServiceProvider } from "ketcher-standalone";
import "ketcher-react/dist/index.css";

// Extend Window interface for TypeScript
declare global {
  interface Window {
    ketcher: any;
  }
}

// Match the example's pattern for getting struct service provider
async function getStructServiceProvider(): Promise<StructServiceProvider> {
  return new StandaloneStructServiceProvider();
}

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
  const [structServiceProvider, setStructServiceProvider] =
    useState<StructServiceProvider | null>(null);

  // Load struct service provider asynchronously (matches official example pattern)
  useEffect(() => {
    getStructServiceProvider().then(setStructServiceProvider);
  }, []);

  // Sync local state with incoming data prop (e.g., after form reset)
  useEffect(() => {
    setCurrentSmiles(data ?? "");
  }, [data]);


  const handleOpenDialog = () => {
    setDialogOpen(true);
  };

  const handleCloseDialog = () => {
    setDialogOpen(false);
  };

  const handleApply = async () => {
    try {
      if (window.ketcher) {
        const smiles = await window.ketcher.getSmiles();
        setCurrentSmiles(smiles);
        handleChange(path, smiles);
      }
    } catch (error) {
      console.error("Error getting SMILES from Ketcher:", error);
    }
    setDialogOpen(false);
  };

  const handleTextChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newSmiles = event.target.value;
    setCurrentSmiles(newSmiles);
    handleChange(path, newSmiles);
  };


  return (
    <Box sx={{ marginBottom: 2 }}>
      <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
        <FormLabel required={required}>{label}</FormLabel>
        {schema.description && (
          <Tooltip title={schema.description} placement="top">
            <IconButton size="small" sx={{ padding: 0.5 }}>
              <HelpOutlineIcon fontSize="small" />
            </IconButton>
          </Tooltip>
        )}
      </Box>
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

      <Dialog
        open={dialogOpen}
        onClose={handleCloseDialog}
        maxWidth="lg"
        fullWidth
      >
        <DialogTitle>Molecular Structure Editor</DialogTitle>
        <DialogContent>
          {!structServiceProvider ? (
            <Box sx={{ height: "400px", display: "flex", alignItems: "center", justifyContent: "center" }}>
              Loading...
            </Box>
          ) : (
              <Box sx={{ height: "400px", width: "100%" }}>
                <Editor
                  errorHandler={(message: string) => {
                    console.error("Ketcher Error:", message);
                  }}
                  staticResourcesUrl={""}
                  onInit={(ketcher) => {
                    window.ketcher = ketcher;
                    ketcher.setMolecule(currentSmiles || "");
                  }}
                  structServiceProvider={structServiceProvider}
                  disableMacromoleculesEditor
                />
              </Box>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseDialog}>Cancel</Button>
          <Button onClick={handleApply} variant="contained" color="primary">
            Apply
          </Button>
        </DialogActions>
      </Dialog>
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
