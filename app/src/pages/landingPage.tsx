import Box from "@mui/material/Box";
import AppBar from "@mui/material/AppBar";
import CssBaseline from "@mui/material/CssBaseline";
import Toolbar from "@mui/material/Toolbar";
import Typography from "@mui/material/Typography";
import IconButton from "@mui/material/IconButton";
import Tooltip from "@mui/material/Tooltip";
import ChatIcon from "@mui/icons-material/Chat";
import CodeIcon from "@mui/icons-material/Code";
import LightModeIcon from "@mui/icons-material/LightMode";
import DarkModeIcon from "@mui/icons-material/DarkMode";
import BrushIcon from "@mui/icons-material/Brush";
import FrameProgressBar from "../components/ProgressBar";
import SideBar from "../components/SideBar";
import RoomManagementMenu from "../components/RoomManagementMenu";

import { useSocketManager } from "../hooks/useSocketManager";
import { useKeyboardShortcuts } from "../hooks/useKeyboardShortcuts";
import MyScene from "../components/Canvas";
import ChatWindow from "../components/ChatWindow";
import ConnectionDialog from "../components/ConnectionDialog";
import { useAppStore } from "../store";
import { useRestJoinManager } from "../hooks/useRestManager";
import { useState } from "react";
import { useColorScheme } from "@mui/material/styles";
import WindowManager from "../components/WindowManager";
import AddPlotButton from "../components/AddPlotButton";
import { useQueryClient } from "@tanstack/react-query";

export default function MainPage() {
  useSocketManager();
  useKeyboardShortcuts();
  useRestJoinManager();

  const { chatOpen, setChatOpen, isDrawing, toggleDrawingMode } = useAppStore();
  const [connectionDialogOpen, setConnectionDialogOpen] = useState(false);
  const { mode, setMode } = useColorScheme();
  const queryClient = useQueryClient();

  const handleToggleColorMode = () => {
    setMode(mode === "light" ? "dark" : "light");
  };

  return (
    <>
      <Box
        sx={{
          display: "flex",
          flexDirection: "column",
          height: "100vh",
          width: "100vw",
          overflow: "hidden",
        }}
      >
        <CssBaseline />

        {/* Header / AppBar */}
        <AppBar
          position="static"
          sx={{ zIndex: (theme) => theme.zIndex.drawer + 1 }}
        >
          <Toolbar>
            <Typography
              variant="h6"
              noWrap
              component="div"
              sx={{ flexGrow: 1 }}
            >
              ZnDraw
            </Typography>
            <Tooltip title={isDrawing ? "Disable drawing mode" : "Enable drawing mode"}>
              <IconButton
                color="inherit"
                aria-label="toggle drawing mode"
                onClick={() => toggleDrawingMode(queryClient)}
                sx={{
                  backgroundColor: isDrawing ? 'rgba(255, 255, 255, 0.2)' : 'transparent',
                }}
              >
                <BrushIcon />
              </IconButton>
            </Tooltip>
            <Tooltip
              title={
                mode === "light"
                  ? "Switch to dark mode"
                  : "Switch to light mode"
              }
            >
              <IconButton
                color="inherit"
                aria-label="toggle color mode"
                onClick={handleToggleColorMode}
              >
                {mode === "light" ? <DarkModeIcon /> : <LightModeIcon />}
              </IconButton>
            </Tooltip>
            <Tooltip title={"Python code connection info"}>
              <IconButton
                color="inherit"
                aria-label="show connection info"
                onClick={() => setConnectionDialogOpen(true)}
              >
                <CodeIcon />
              </IconButton>
            </Tooltip>
            <Tooltip title={"Toggle chat window"}>
              <IconButton
                color="inherit"
                aria-label="toggle chat"
                onClick={() => setChatOpen(!chatOpen)}
              >
                <ChatIcon />
              </IconButton>
            </Tooltip>
            <AddPlotButton />
            <RoomManagementMenu />
          </Toolbar>
        </AppBar>

        {/* Main content row with sidebar and center area */}
        <Box sx={{ display: "flex", flexGrow: 1, minHeight: 0 }}>
          <SideBar />

          {/* Main Content Area with drag boundary */}
          <Box
            component="main"
            sx={{
              flexGrow: 1,
              display: "flex",
              flexDirection: "column",
              minWidth: 0,
            }}
          >
            {/* THIS IS THE CRUCIAL DRAG BOUNDARY CONTAINER */}
            <Box
              className="drag-boundary-container"
              sx={{
                flexGrow: 1,
                position: "relative",
                overflow: "hidden",
                display: "flex",
                flexDirection: "column",
              }}
            >
              <MyScene />
              <WindowManager />
            </Box>
          </Box>
        </Box>

        {/* Progress bar spans full width outside the drag boundary */}
        <FrameProgressBar />

        {/* Chat Window */}
        <ChatWindow open={chatOpen} onClose={() => setChatOpen(false)} />
      </Box>

      <ConnectionDialog
        open={connectionDialogOpen}
        onClose={() => setConnectionDialogOpen(false)}
      />
    </>
  );
}
