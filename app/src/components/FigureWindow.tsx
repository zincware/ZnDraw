import { Rnd } from 'react-rnd';
import Plot from 'react-plotly.js';
import { useWindowManagerStore } from '../stores/windowManagerStore';
import { useFigure } from '../hooks/useFigures';

// MUI Imports
import { Box, Paper, Typography, IconButton, CircularProgress } from '@mui/material';
import { Close as CloseIcon, ErrorOutline as ErrorIcon } from '@mui/icons-material';

interface FigureWindowProps {
  figureKey: string;
}

function FigureWindow({ figureKey }: FigureWindowProps) {
  const windowState = useWindowManagerStore((state) => state.openWindows[figureKey]);
  const { updateWindowState, closeWindow, bringToFront } = useWindowManagerStore();
  const { data: figureResponse, isLoading, isError, error } = useFigure(figureKey);

  if (!windowState) return null;

  // Function to render the content based on the fetch state
  const renderContent = () => {
    if (isLoading) {
      return (
        <Box display="flex" justifyContent="center" alignItems="center" height="100%">
          <CircularProgress />
        </Box>
      );
    }

    if (isError) {
      return (
        <Box
          display="flex"
          flexDirection="column"
          justifyContent="center"
          alignItems="center"
          height="100%"
          color="error.main"
        >
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
                // Your backend stores the plotly figure as a JSON string
                const plotlyJson = JSON.parse(figureData.data);
                return (
                    <Plot
                        data={plotlyJson.data}
                        layout={{ ...plotlyJson.layout, autosize: true }}
                        useResizeHandler={true}
                        style={{ width: '100%', height: '100%' }}
                    />
                );
            } catch (e) {
                return <Typography color="error">Invalid Plotly JSON format.</Typography>
            }
        }
        // Fallback for other figure types like 'line' or 'circle'
        return <pre>{JSON.stringify(figureData, null, 2)}</pre>;
    }

    return null;
  };

  return (
    <Rnd
      size={{ width: windowState.width, height: windowState.height }}
      position={{ x: windowState.x, y: windowState.y }}
      minWidth={300}
      minHeight={250}
      style={{ zIndex: windowState.zIndex }}
      dragHandleClassName="window-header" // This class name is used by Rnd to identify the drag handle
      onDragStart={() => bringToFront(figureKey)}
      onDragStop={(_, d) => {
        updateWindowState(figureKey, { x: d.x, y: d.y });
      }}
      onResizeStop={(_, __, ref, ___, position) => {
        updateWindowState(figureKey, {
          width: parseInt(ref.style.width, 10),
          height: parseInt(ref.style.height, 10),
          ...position,
        });
      }}
      bounds="parent"
    >
      <Paper 
        elevation={8}
        sx={{
            width: '100%',
            height: '100%',
            display: 'flex',
            flexDirection: 'column',
            overflow: 'hidden', // Ensures content doesn't spill out
        }}
      >
        {/* Header / Title Bar */}
        <Box
          className="window-header" // Critical for react-rnd dragging
          onMouseDown={() => bringToFront(figureKey)}
          sx={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            bgcolor: 'primary.main',
            color: 'primary.contrastText',
            py: 0.5,
            px: 2,
            cursor: 'grab',
          }}
        >
          <Typography variant="subtitle2" component="h2" noWrap sx={{ flexGrow: 1 }}>
            Figure: {figureKey}
          </Typography>
          <IconButton onClick={() => closeWindow(figureKey)} size="small" color="inherit">
            <CloseIcon fontSize="small" />
          </IconButton>
        </Box>

        {/* Content Area */}
        <Box
          sx={{
            flexGrow: 1,
            p: 1, // Add some padding around the content
            overflow: 'auto', // Allow content to scroll if it overflows
            position: 'relative' // For centering loading/error states
          }}
        >
          {renderContent()}
        </Box>
      </Paper>
    </Rnd>
  );
}

export default FigureWindow;