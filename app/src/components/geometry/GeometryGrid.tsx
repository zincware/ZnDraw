import { useMemo, useState } from "react";
import {
  DataGrid,
  GridColDef,
  GridActionsCellItem,
  GridRowParams,
} from "@mui/x-data-grid";
import { Box, TextField, InputAdornment } from "@mui/material";
import EditIcon from "@mui/icons-material/Edit";
import DeleteIcon from "@mui/icons-material/Delete";
import SearchIcon from "@mui/icons-material/Search";
import { useGeometryStore } from "../../stores/geometryStore";
import DeleteConfirmDialog from "./DeleteConfirmDialog";
import { useDeleteGeometry } from "../../hooks/useGeometries";
import { useAppStore } from "../../store";

interface GeometryGridProps {
  geometries: Array<{ key: string; type: string }>;
}

const GeometryGrid = ({ geometries }: GeometryGridProps) => {
  const { roomId } = useAppStore();
  const { setMode, setSelectedKey, searchFilter, setSearchFilter } =
    useGeometryStore();
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [geometryToDelete, setGeometryToDelete] = useState<string | null>(null);

  const { mutate: deleteGeometry } = useDeleteGeometry();

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
