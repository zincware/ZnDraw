// src/components/jsonforms-renderers/SmilesEditDialog.tsx
// Shared Ketcher-based SMILES editor dialog

import { useState, useEffect } from "react";
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  Box,
} from "@mui/material";
import { Editor } from "ketcher-react";
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

interface SmilesEditDialogProps {
  open: boolean;
  initialSmiles: string;
  title?: string;
  onSave: (smiles: string) => void;
  onCancel: () => void;
}

/**
 * Reusable Ketcher-based SMILES editor dialog.
 * Can be used by any component that needs molecular structure editing.
 */
export default function SmilesEditDialog({
  open,
  initialSmiles,
  title = "Molecular Structure Editor",
  onSave,
  onCancel,
}: SmilesEditDialogProps) {
  const [structServiceProvider, setStructServiceProvider] =
    useState<StructServiceProvider | null>(null);

  // Load struct service provider asynchronously
  useEffect(() => {
    getStructServiceProvider().then(setStructServiceProvider);
  }, []);

  const handleApply = async () => {
    try {
      if (window.ketcher) {
        const smiles = await window.ketcher.getSmiles();
        onSave(smiles);
      }
    } catch (error) {
      console.error("Error getting SMILES from Ketcher:", error);
    }
  };

  return (
    <Dialog
      open={open}
      onClose={onCancel}
      maxWidth="lg"
      fullWidth
    >
      <DialogTitle>{title}</DialogTitle>
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
                ketcher.setMolecule(initialSmiles || "");
              }}
              structServiceProvider={structServiceProvider}
              disableMacromoleculesEditor
            />
          </Box>
        )}
      </DialogContent>
      <DialogActions>
        <Button onClick={onCancel}>Cancel</Button>
        <Button onClick={handleApply} variant="contained" color="primary">
          Apply
        </Button>
      </DialogActions>
    </Dialog>
  );
}
