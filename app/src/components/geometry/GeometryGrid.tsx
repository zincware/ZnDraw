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
  const { roomId, geometries: geometriesData, activeCurveForDrawing, setActiveCurveForDrawing } = useAppStore();
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
      clientId: null,
      key,
      geometryType: geometry.type,
      geometryData: {
        active: !currentActive,
      },
    });
  };

  const handleDrawingTargetToggle = (key: string, isCurve: boolean) => {
    if (!isCurve) return;

    // Toggle: if already selected, deselect; otherwise select
    if (activeCurveForDrawing === key) {
      setActiveCurveForDrawing(null);
    } else {
      setActiveCurveForDrawing(key);
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
      field: "drawingTarget",
      headerName: "Draw",
      width: 70,
      renderCell: (params: GridRenderCellParams) => {
        const isCurve = params.row.type === "Curve";
        const isSelected = activeCurveForDrawing === params.row.key;

        if (!isCurve) {
          return null;
        }

        return (
          <Tooltip title={isSelected ? "Stop drawing on this curve" : "Draw on this curve"}>
            <IconButton
              size="small"
              onClick={(e) => {
                e.stopPropagation();
                handleDrawingTargetToggle(params.row.key, isCurve);
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
      headerName: "Active",
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
