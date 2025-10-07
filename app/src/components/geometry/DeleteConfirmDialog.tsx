import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogContentText,
  DialogActions,
  Button,
} from "@mui/material";

interface DeleteConfirmDialogProps {
  open: boolean;
  geometryKey: string | null;
  onConfirm: () => void;
  onCancel: () => void;
}

const DeleteConfirmDialog = ({
  open,
  geometryKey,
  onConfirm,
  onCancel,
}: DeleteConfirmDialogProps) => {
  return (
    <Dialog open={open} onClose={onCancel}>
      <DialogTitle>Delete Geometry</DialogTitle>
      <DialogContent>
        <DialogContentText>
          Are you sure you want to delete the geometry "{geometryKey}"? This
          action cannot be undone.
        </DialogContentText>
      </DialogContent>
      <DialogActions>
        <Button onClick={onCancel}>Cancel</Button>
        <Button onClick={onConfirm} color="error" variant="contained">
          Delete
        </Button>
      </DialogActions>
    </Dialog>
  );
};

export default DeleteConfirmDialog;
