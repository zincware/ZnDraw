import { useMemo } from "react";
import {
  Box,
  Typography,
  Divider,
  Button,
  CircularProgress,
} from "@mui/material";
import AddIcon from "@mui/icons-material/Add";
import { useGeometryStore } from "../../stores/geometryStore";
import { useGeometriesList } from "../../hooks/useGeometries";
import { useAppStore } from "../../store";
import GeometryGrid from "./GeometryGrid";
import GeometryForm from "./GeometryForm";

const GeometryPanel = () => {
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const { mode, setMode, resetForm } = useGeometryStore();

  const {
    data: geometriesData,
    isLoading,
    isError,
  } = useGeometriesList(roomId);

  const geometries = useMemo(() => {
    if (!geometriesData?.geometries || typeof geometriesData.geometries !== 'object') {
      return [];
    }

    // geometriesData.geometries is an object with geometry keys as properties
    // Each property has {type: string, data: object} structure
    return Object.keys(geometriesData.geometries).map((key: string) => ({
      key,
      type: geometriesData.geometries[key]?.type || "unknown",
    }));
  }, [geometriesData]);

  const handleCreate = () => {
    resetForm();
    setMode("create");
  };

  if (!roomId) {
    return <Typography sx={{ p: 2 }}>Joining room...</Typography>;
  }

  if (isLoading) {
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

  if (isError) {
    return (
      <Typography color="error" sx={{ p: 2 }}>
        Failed to load geometries.
      </Typography>
    );
  }

  return (
    <Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
      {mode === "list" && (
        <>
          <Typography variant="h6" sx={{ p: 2, pb: 1, flexShrink: 0 }}>
            Geometries
          </Typography>
          <Divider />
          <Box sx={{ p: 2 }}>
            <Button
              variant="contained"
              startIcon={<AddIcon />}
              onClick={handleCreate}
              fullWidth
            >
              Add Geometry
            </Button>
          </Box>
          <Box sx={{ flexGrow: 1, minHeight: 0 }}>
            <GeometryGrid geometries={geometries} />
          </Box>
        </>
      )}
      {(mode === "create" || mode === "edit") && <GeometryForm />}
    </Box>
  );
};

export default GeometryPanel;
