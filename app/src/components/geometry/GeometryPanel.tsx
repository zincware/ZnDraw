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
  const { roomId } = useAppStore();
  const { mode, setMode, resetForm } = useGeometryStore();

  const {
    data: geometriesData,
    isLoading,
    isError,
  } = useGeometriesList(roomId);

  const geometries = useMemo(() => {
    if (!geometriesData?.geometries || !Array.isArray(geometriesData.geometries)) {
      return [];
    }

    // Get geometries from the store (which has type info)
    const geometriesObj = useAppStore.getState().geometries;

    return geometriesData.geometries.map((key: string) => ({
      key,
      type: geometriesObj[key]?.type || "unknown",
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
