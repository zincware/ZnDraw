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
    Select, 
    MenuItem, 
    FormControl 
} from '@mui/material';
import { SelectChangeEvent } from '@mui/material/Select';
import { Close as CloseIcon, ErrorOutline as ErrorIcon } from '@mui/icons-material';

interface FigureWindowProps {
  windowId: string; // The component is now identified by its own stable ID
}

function FigureWindow({ windowId }: FigureWindowProps) {
  // --- STATE & DATA ---
  
  // 1. Get the window's own state (position, size, and current figureKey) from the store
  const windowInstance = useWindowManagerStore((state) => state.openWindows[windowId]);
  
  // 2. Get the actions from the store to manipulate this and other windows
  const { 
    updateWindowState, 
    closeWindow, 
    bringToFront, 
    changeFigureInWindow 
  } = useWindowManagerStore();
  
  // 3. Fetch the actual plot data for the figure this window is currently displaying
  // The 'enabled' flag prevents fetching if the window instance isn't ready
  const { 
    data: figureResponse, 
    isLoading, 
    isError, 
    error 
  } = useFigure(windowInstance?.figureKey, { enabled: !!windowInstance });
  
  // 4. Fetch the list of ALL available plots to populate the dropdown menu
  const { data: allFiguresResponse } = useFigureList();

  // --- GUARDS ---
  
  // If the window was closed, its instance will be removed from the store, so we render nothing.
  if (!windowInstance) {
    return null;
  }

  // --- HANDLERS ---
  
  const handleFigureChange = (event: SelectChangeEvent<string>) => {
    const newFigureKey = event.target.value;
    // Call the store action to update which figure this window is showing
    changeFigureInWindow(windowId, newFigureKey);
  };
  
  // --- RENDER LOGIC ---

  /**
   * Renders the content of the window based on the data fetching status.
   */
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
          textAlign="center"
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
                // The backend stores the plotly figure as a JSON string, so we parse it
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
                return <Typography color="error">Invalid Plotly JSON format.</Typography>;
            }
        }
        // Fallback for other figure types like 'line' or 'circle'
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
      dragHandleClassName="window-header" // This class name is used by Rnd to identify the drag handle
      onDragStart={() => bringToFront(windowId)}
      onDragStop={(_, d) => {
        updateWindowState(windowId, { x: d.x, y: d.y });
      }}
      onResizeStop={(_, __, ref, ___, position) => {
        updateWindowState(windowId, {
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
          <Typography variant="subtitle2" component="h2" noWrap>Figure:</Typography>
          
          <FormControl variant="standard" size="small" sx={{ flexGrow: 1, minWidth: 120 }}>
            <Select
              value={windowInstance.figureKey}
              onChange={handleFigureChange}
              // Styling to make the select blend with the dark header
              sx={{
                color: 'inherit',
                '& .MuiSelect-icon': { color: 'inherit' },
                '&:before': { borderBottom: '1px solid rgba(255, 255, 255, 0.7)' },
                '&:hover:not(.Mui-disabled):before': { borderBottom: '1px solid #fff' },
              }}
            >
              {allFiguresResponse?.figures.map((key) => (
                <MenuItem key={key} value={key}>
                  {key}
                </MenuItem>
              ))}
            </Select>
          </FormControl>
          
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
            position: 'relative', // For centering loading/error states
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