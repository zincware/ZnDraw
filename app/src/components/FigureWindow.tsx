import { Rnd } from 'react-rnd';
import Plot from 'react-plotly.js';
import { useWindowManagerStore } from '../stores/windowManagerStore';
import { useFigure, useFigureList } from '../hooks/useFigures';

// MUI Imports
import { 
    Box, 
    Paper, 
    Typography, 
    IconButton, 
    CircularProgress, 
    TextField,
    Autocomplete
} from '@mui/material';
import { Close as CloseIcon, ErrorOutline as ErrorIcon } from '@mui/icons-material';

interface FigureWindowProps {
  windowId: string;
}

function FigureWindow({ windowId }: FigureWindowProps) {
  const windowInstance = useWindowManagerStore((state) => state.openWindows[windowId]);
  const { 
    updateWindowState, 
    closeWindow, 
    bringToFront, 
    changeFigureInWindow 
  } = useWindowManagerStore();
  const { 
    data: figureResponse, 
    isLoading, 
    isError, 
    error 
  } = useFigure(windowInstance?.figureKey, { enabled: !!windowInstance });
  const { data: allFiguresResponse } = useFigureList();

  if (!windowInstance) {
    return null;
  }

  const handleAutocompleteChange = (event: React.SyntheticEvent, newValue: string | null) => {
    if (newValue) {
      changeFigureInWindow(windowId, newValue);
    }
  };
  
  const renderContent = () => {
    if (isLoading) {
      return <Box display="flex" justifyContent="center" alignItems="center" height="100%"><CircularProgress /></Box>;
    }
    if (isError) {
      return (
        <Box display="flex" flexDirection="column" justifyContent="center" alignItems="center" height="100%" color="error.main" textAlign="center">
          <ErrorIcon sx={{ fontSize: 40, mb: 1 }} />
          <Typography variant="body1">Error loading figure.</Typography>
          <Typography variant="caption">{error?.message}</Typography>
        </Box>
      );
    }
    if (figureResponse?.figure) {
        const figureData = figureResponse.figure;
        if (figureData.type === 'plotly' && figureData.data) {
            try {
                const plotlyJson = JSON.parse(figureData.data);
                return <Plot data={plotlyJson.data} layout={{ ...plotlyJson.layout, autosize: true }} useResizeHandler={true} style={{ width: '100%', height: '100%' }} />;
            } catch (e) {
                return <Typography color="error">Invalid Plotly JSON format.</Typography>;
            }
        }
        return <pre style={{ margin: 0 }}>{JSON.stringify(figureData, null, 2)}</pre>;
    }
    return null;
  };

  return (
    <Rnd
      size={{ width: windowInstance.width, height: windowInstance.height }}
      position={{ x: windowInstance.x, y: windowInstance.y }}
      minWidth={350}
      minHeight={300}
      style={{ zIndex: windowInstance.zIndex }}
      dragHandleClassName="window-header"
      onDragStart={() => bringToFront(windowId)}
      onDragStop={(_, d) => { updateWindowState(windowId, { x: d.x, y: d.y }); }}
      onResizeStop={(_, __, ref, ___, position) => {
        updateWindowState(windowId, {
          width: parseInt(ref.style.width, 10),
          height: parseInt(ref.style.height, 10),
          ...position,
        });
      }}
      bounds=".drag-boundary-container"
    >
      <Paper 
        elevation={8}
        sx={{
            width: '100%',
            height: '100%',
            display: 'flex',
            flexDirection: 'column',
            overflow: 'hidden',
            border: '1px solid rgba(0,0,0,0.2)'
        }}
      >
        {/* Header / Title Bar */}
        <Box
          className="window-header"
          onMouseDown={() => bringToFront(windowId)}
          sx={{
            display: 'flex',
            alignItems: 'center',
            gap: 2,
            bgcolor: 'primary.main',
            color: 'primary.contrastText',
            py: 0.5,
            px: 2,
            cursor: 'grab',
          }}
        >
          {/* --- FIXES APPLIED HERE --- */}
          <Autocomplete
            options={allFiguresResponse?.figures || []}
            value={windowInstance.figureKey}
            onChange={handleAutocompleteChange}
            fullWidth
            disableClearable
            size="small"
            // FIX 1: This is the most important change. It stops the click from
            // starting a drag action on the window header.
            onMouseDown={(event) => event.stopPropagation()}
            // FIX 2: This renders the dropdown list within the Rnd component,
            // preventing z-index and positioning issues.
            disablePortal
            renderInput={(params) => (
              <TextField
                {...params}
                label="Figure"
                variant="standard"
                // Simplified styling for better readability and performance
                sx={{
                  '& .MuiInput-underline:before': { borderBottomColor: 'rgba(255, 255, 255, 0.7)' },
                  '& .MuiInputBase-input': { color: 'primary.contrastText' },
                  '& .MuiInputLabel-root': { color: 'rgba(255, 255, 255, 0.7)' },
                  '& .MuiSvgIcon-root': { color: 'primary.contrastText' }
                }}
              />
            )}
          />
          
          <IconButton onClick={() => closeWindow(windowId)} size="small" color="inherit">
            <CloseIcon fontSize="small" />
          </IconButton>
        </Box>

        {/* Content Area */}
        <Box
          sx={{
            flexGrow: 1,
            p: 1,
            overflow: 'auto',
            position: 'relative',
            backgroundColor: '#fff'
          }}
        >
          {renderContent()}
        </Box>
      </Paper>
    </Rnd>
  );
}

export default FigureWindow;