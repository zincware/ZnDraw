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
  CircularProgress,
  Alert,
} from "@mui/material";
import { rankWith, schemaMatches, type ControlProps } from "@jsonforms/core";
import { useState, useEffect, useCallback } from "react";
import { Editor } from "ketcher-react";
import { StructServiceProvider } from "ketcher-core";
import { StandaloneStructServiceProvider } from "ketcher-standalone";
import { throttle } from "lodash";
import { convertMoleculeToImage } from "../../myapi/client";
import { LAYOUT_CONSTANTS } from "../../constants/layout";
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
  const [moleculeImage, setMoleculeImage] = useState<string | null>(null);
  const [imageLoading, setImageLoading] = useState(false);
  const [imageError, setImageError] = useState<string | null>(null);

  // Load struct service provider asynchronously (matches official example pattern)
  useEffect(() => {
    getStructServiceProvider().then(setStructServiceProvider);
  }, []);

  // Sync local state with incoming data prop (e.g., after form reset)
  useEffect(() => {
    setCurrentSmiles(data ?? "");
  }, [data]);

  // Throttled function to fetch molecule image
  const fetchMoleculeImage = useCallback(
    throttle(async (smiles: string) => {
      if (!smiles || smiles.trim() === "") {
        setMoleculeImage(null);
        setImageError(null);
        return;
      }

      setImageLoading(true);
      setImageError(null);

      try {
        const response = await convertMoleculeToImage({
          type: "smiles",
          data: smiles,
        });
        setMoleculeImage(response.image);
      } catch (error: any) {
        console.error("Error fetching molecule image:", error);
        setImageError(error.response?.data?.error || "Failed to generate image");
        setMoleculeImage(null);
      } finally {
        setImageLoading(false);
      }
    }, 500), // Throttle to 500ms
    []
  );

  // Fetch image when currentSmiles changes
  useEffect(() => {
    fetchMoleculeImage(currentSmiles);
  }, [currentSmiles, fetchMoleculeImage]);


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

      {/* Molecule Image Preview */}
      {currentSmiles && (
        <Box sx={{ marginTop: 2, display: "flex", justifyContent: "center", alignItems: "center", minHeight: 100 }}>
          {imageLoading && <CircularProgress size={40} />}
          {imageError && (
            <Alert severity="error" sx={{ width: "100%" }}>
              {imageError}
            </Alert>
          )}
          {moleculeImage && !imageLoading && !imageError && (
            <Box
              component="img"
              src={moleculeImage}
              alt="Molecule structure"
              sx={{
                width: "100%",
                maxWidth: `calc(100vw - ${LAYOUT_CONSTANTS.PRIMARY_DRAWER_WIDTH}px - 48px)`,
                maxHeight: "300px",
                objectFit: "contain",
                border: "1px solid #e0e0e0",
                borderRadius: 1,
                padding: 1,
                backgroundColor: "white",
              }}
            />
          )}
        </Box>
      )}

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
