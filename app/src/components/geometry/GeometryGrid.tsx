import { useMemo, useState } from "react";
import {
  DataGrid,
  GridColDef,
  GridActionsCellItem,
  GridRowParams,
  GridRenderCellParams,
} from "@mui/x-data-grid";
import { Box, TextField, InputAdornment, Switch, IconButton, Tooltip } from "@mui/material";
import EditIcon from "@mui/icons-material/Edit";
import DeleteIcon from "@mui/icons-material/Delete";
import SearchIcon from "@mui/icons-material/Search";
import BrushIcon from "@mui/icons-material/Brush";
import RadioButtonUncheckedIcon from "@mui/icons-material/RadioButtonUnchecked";
import RadioButtonCheckedIcon from "@mui/icons-material/RadioButtonChecked";
import { useGeometryStore } from "../../stores/geometryStore";
import DeleteConfirmDialog from "./DeleteConfirmDialog";
import { useDeleteGeometry, useCreateGeometry } from "../../hooks/useGeometries";
import { useAppStore } from "../../store";

interface GeometryGridProps {
  geometries: Array<{ key: string; type: string }>;
}

const GeometryGrid = ({ geometries }: GeometryGridProps) => {
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const geometriesData = useAppStore((state) => state.geometries);
  const activeCurveForDrawing = useAppStore((state) => state.activeCurveForDrawing);
  const setActiveCurveForDrawing = useAppStore((state) => state.setActiveCurveForDrawing);
  const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);
  const attachToCamera = useAppStore((state) => state.attachToCamera);
  const detachFromCamera = useAppStore((state) => state.detachFromCamera);
  const { setMode, setSelectedKey, searchFilter, setSearchFilter } =
    useGeometryStore();
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [geometryToDelete, setGeometryToDelete] = useState<string | null>(null);

  const { mutate: deleteGeometry } = useDeleteGeometry();
  const { mutate: updateGeometry } = useCreateGeometry();

  const handleEdit = (key: string) => {
    setSelectedKey(key);
    setMode("edit");
  };

  const handleDeleteClick = (key: string) => {
    setGeometryToDelete(key);
    setDeleteDialogOpen(true);
  };

  const handleDeleteConfirm = () => {
    if (geometryToDelete && roomId) {
      deleteGeometry({ roomId, key: geometryToDelete });
    }
    setDeleteDialogOpen(false);
    setGeometryToDelete(null);
  };

  const handleDeleteCancel = () => {
    setDeleteDialogOpen(false);
    setGeometryToDelete(null);
  };

  const handleActiveToggle = (key: string, currentActive: boolean) => {
    const geometry = geometriesData[key];
    if (!geometry || !roomId) return;

    // Update the geometry with the new active state
    updateGeometry({
      roomId,
      userName: null,
      key,
      geometryType: geometry.type,
      geometryData: {
        active: !currentActive,
      },
    });
  };

  const handleSelectionToggle = (key: string, geometryType: string) => {
    if (geometryType === "Curve") {
      // Toggle curve drawing target
      if (activeCurveForDrawing === key) {
        setActiveCurveForDrawing(null);
      } else {
        setActiveCurveForDrawing(key);
      }
    } else if (geometryType === "Camera") {
      // Toggle camera attachment
      if (attachedCameraKey === key) {
        detachFromCamera();
      } else {
        attachToCamera(key);
      }
    }
  };

  const columns: GridColDef[] = [
    {
      field: "key",
      headerName: "Key",
      flex: 1,
      minWidth: 150,
    },
    {
      field: "type",
      headerName: "Type",
      flex: 1,
      minWidth: 120,
    },
    {
      field: "selection",
      headerName: "Active",
      width: 70,
      renderCell: (params: GridRenderCellParams) => {
        const geometryType = params.row.type;
        const isCurve = geometryType === "Curve";
        const isCamera = geometryType === "Camera";

        // Only show for Curves and Cameras
        if (!isCurve && !isCamera) {
          return null;
        }

        const isSelected = isCurve
          ? activeCurveForDrawing === params.row.key
          : attachedCameraKey === params.row.key;

        const tooltipText = isCurve
          ? (isSelected ? "Deselect this curve for drawing" : "Select this curve for drawing")
          : (isSelected ? "Detach from this camera" : "Attach to this camera");

        return (
          <Tooltip title={tooltipText}>
            <IconButton
              size="small"
              onClick={(e) => {
                e.stopPropagation();
                handleSelectionToggle(params.row.key, geometryType);
              }}
              color={isSelected ? "primary" : "default"}
            >
              {isSelected ? <RadioButtonCheckedIcon /> : <RadioButtonUncheckedIcon />}
            </IconButton>
          </Tooltip>
        );
      },
    },
    {
      field: "active",
      headerName: "Visible",
      width: 80,
      renderCell: (params: GridRenderCellParams) => {
        const isActive = params.row.active !== false;
        return (
          <Switch
            checked={isActive}
            onChange={(e) => {
              e.stopPropagation();
              handleActiveToggle(params.row.key, isActive);
            }}
            onClick={(e) => {
              e.stopPropagation();
            }}
            size="small"
          />
        );
      },
    },
    {
      field: "actions",
      type: "actions",
      headerName: "Actions",
      width: 100,
      getActions: (params: GridRowParams) => [
        <GridActionsCellItem
          icon={<EditIcon />}
          label="Edit"
          onClick={() => handleEdit(params.row.key)}
          showInMenu={false}
        />,
        <GridActionsCellItem
          icon={<DeleteIcon />}
          label="Delete"
          onClick={() => handleDeleteClick(params.row.key)}
          showInMenu={false}
        />,
      ],
    },
  ];

  const filteredGeometries = useMemo(() => {
    if (!searchFilter) return geometries;
    const lowerFilter = searchFilter.toLowerCase();
    return geometries.filter(
      (g) =>
        g.key.toLowerCase().includes(lowerFilter) ||
        g.type.toLowerCase().includes(lowerFilter)
    );
  }, [geometries, searchFilter]);

  const rows = filteredGeometries.map((g, index) => ({
    id: index,
    key: g.key,
    type: g.type,
    active: geometriesData[g.key]?.data?.active !== false, // Default to true if undefined
  }));

  return (
    <Box sx={{ height: "100%", display: "flex", flexDirection: "column" }}>
      <Box sx={{ p: 2, pb: 1 }}>
        <TextField
          fullWidth
          size="small"
          placeholder="Search geometries..."
          value={searchFilter}
          onChange={(e) => setSearchFilter(e.target.value)}
          slotProps={{
            input: {
              startAdornment: (
                <InputAdornment position="start">
                  <SearchIcon />
                </InputAdornment>
              ),
            },
          }}
        />
      </Box>
      <Box sx={{ flexGrow: 1, px: 2, pb: 2 }}>
        <DataGrid
          rows={rows}
          columns={columns}
          density="compact"
          disableRowSelectionOnClick
          sx={{
            border: "1px solid rgba(0, 0, 0, 0.12)",
            "& .MuiDataGrid-cell:focus": {
              outline: "none",
            },
            "& .MuiDataGrid-row:hover": {
              cursor: "pointer",
            },
          }}
          onRowClick={(params) => handleEdit(params.row.key)}
        />
      </Box>
      <DeleteConfirmDialog
        open={deleteDialogOpen}
        geometryKey={geometryToDelete}
        onConfirm={handleDeleteConfirm}
        onCancel={handleDeleteCancel}
      />
    </Box>
  );
};

export default GeometryGrid;
