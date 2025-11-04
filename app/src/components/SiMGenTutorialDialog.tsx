import {
  Dialog,
  DialogTitle,
  DialogContent,
  IconButton,
} from '@mui/material';
import CloseIcon from '@mui/icons-material/Close';

interface SiMGenTutorialDialogProps {
  open: boolean;
  onClose: () => void;
  url: string;
}

export default function SiMGenTutorialDialog({ open, onClose, url }: SiMGenTutorialDialogProps) {
  return (
    <Dialog
      open={open}
      onClose={onClose}
      maxWidth="lg"
      fullWidth
      slotProps={{
        paper: {
          sx: {
            height: '90vh',
          },
        },
      }}
    >
      <DialogTitle>
        SiMGen Tutorial
        <IconButton
          aria-label="close"
          onClick={onClose}
          sx={{
            position: 'absolute',
            right: 8,
            top: 8,
            color: (theme) => theme.palette.grey[500],
          }}
        >
          <CloseIcon />
        </IconButton>
      </DialogTitle>
      <DialogContent sx={{ padding: 0 }}>
        <iframe
          src={url}
          title="SiMGen Tutorial"
          style={{
            width: '100%',
            height: '100%',
            border: 'none',
          }}
        />
      </DialogContent>
    </Dialog>
  );
}
