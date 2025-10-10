import { Rnd } from "react-rnd";
import Plot from "react-plotly.js";
import { useWindowManagerStore } from "../stores/windowManagerStore";
import { useFigure, useFigureList } from "../hooks/useFigures";
import { useAppStore } from "../store";
import { useAtomicFrameSet } from "../hooks/useAtomicFrameSet";

// MUI Imports
import {
  Box,
  Paper,
  Typography,
  IconButton,
  CircularProgress,
  TextField,
  Autocomplete,
} from "@mui/material";
import {
  Close as CloseIcon,
  ErrorOutline as ErrorIcon,
} from "@mui/icons-material";

interface FigureWindowProps {
  windowId: string;
}

function FigureWindow({ windowId }: FigureWindowProps) {
  const windowInstance = useWindowManagerStore(
    (state) => state.openWindows[windowId],
  );
  const { updateWindowState, closeWindow, bringToFront, changeFigureInWindow } =
    useWindowManagerStore();

  const {
    data: figureResponse,
    isLoading,
    isError,
    error,
  } = useFigure(windowInstance?.figureKey, { enabled: !!windowInstance });
  const { data: allFiguresResponse } = useFigureList();

  const { setSelection, setFrameSelection } = useAppStore();
  const setFrameAtomic = useAtomicFrameSet();

  if (!windowInstance) {
    return null;
  }

  const onPlotClick = ({ points }: { points: any[] }) => {
    if (!points || points.length === 0) {
      return;
    }
    if (points[0]?.customdata && points[0].customdata[0] != null) {
      setFrameAtomic(points[0].customdata[0]);
    }
    if (points[0]?.customdata && points[0].customdata[1] != null) {
      setSelection([points[0].customdata[1]]);
    }
	};

  const onPlotSelected = (event: any) => {
		if (!event || !event.points) {
			return;
		}
		if (event.points.length === 0) {
			// This is triggered once the plot is re-rendered. We want to keep the selection here.
			return;
		}
		const selectedFrames = event.points.map((point: any) =>
			point.customdata ? point.customdata[0] : point.pointIndex,
		);
		// for all points.customdata[0] == step collect the points.customdata[1] and set selectedIds if customdata[1] is available
		const selectedIds = new Set<number>(
			event.points
				.filter((point: any) => point.customdata?.[1])
				.map((point: any) => point.customdata[1]),
		);
		if (selectedIds.size > 0) {
			setSelection(Array.from(selectedIds));
		}

		setFrameSelection(selectedFrames);
	};

  const onPlotDeselect = () => {
		setFrameSelection([]);
	};

  const handleAutocompleteChange = (
    event: React.SyntheticEvent,
    newValue: string | null,
  ) => {
    if (newValue) {
      changeFigureInWindow(windowId, newValue);
    }
  };

  const renderContent = () => {
    if (isLoading) {
      return (
        <Box
          display="flex"
          justifyContent="center"
          alignItems="center"
          height="100%"
        >
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
      if (figureData.type === "plotly" && figureData.data) {
        try {
          const plotlyJson = JSON.parse(figureData.data);
          return (
            <Plot
              data={plotlyJson.data}
              layout={{ ...plotlyJson.layout, autosize: true }}
              useResizeHandler={true}
              style={{ width: "100%", height: "100%" }}
              onClick={onPlotClick}
              onSelected={onPlotSelected}
              onDeselect={onPlotDeselect}
            />
          );
        } catch (e) {
          return (
            <Typography color="error">Invalid Plotly JSON format.</Typography>
          );
        }
      }
      return (
        <pre style={{ margin: 0 }}>{JSON.stringify(figureData, null, 2)}</pre>
      );
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
      // Use the same class name as the chat for consistency
      dragHandleClassName="drag-handle"
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
      bounds=".drag-boundary-container"
    >
      {/* --- STYLE CHANGE: Elevation and layout match ChatWindow --- */}
      <Paper
        elevation={4}
        sx={{
          width: "100%",
          height: "100%",
          display: "flex",
          flexDirection: "column",
          overflow: "hidden",
        }}
      >
        {/* Header / Title Bar */}
        <Box
          // Use the same class name as the chat for consistency
          className="drag-handle"
          onMouseDown={() => bringToFront(windowId)}
          // --- STYLE CHANGE: Header styles now match ChatWindow ---
          sx={{
            p: 1,
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
            borderBottom: 1,
            borderColor: "divider",
            cursor: "move",
          }}
        >
          {/* We now use a flex layout with space-between like the chat window */}
          <Autocomplete
            options={allFiguresResponse?.figures || []}
            value={windowInstance.figureKey}
            onChange={handleAutocompleteChange}
            disableClearable
            size="small"
            onMouseDown={(event) => event.stopPropagation()}
            disablePortal
            // --- STYLE CHANGE: Autocomplete now has a larger flex footprint ---
            sx={{ flexGrow: 1, minWidth: 150, mr: 2 }}
            renderInput={(params) => (
              <TextField
                {...params}
                label="Figure"
                // --- STYLE CHANGE: Use outlined variant for a cleaner look ---
                variant="outlined"
                size="small"
              />
            )}
          />

          <IconButton onClick={() => closeWindow(windowId)} size="small">
            <CloseIcon />
          </IconButton>
        </Box>

        {/* Content Area */}
        <Box
          sx={{
            flexGrow: 1,
            p: 1,
            overflow: "auto",
            position: "relative",
            backgroundColor: "background.paper", // Use theme background color
          }}
        >
          {renderContent()}
        </Box>
      </Paper>
    </Rnd>
  );
}

export default FigureWindow;
