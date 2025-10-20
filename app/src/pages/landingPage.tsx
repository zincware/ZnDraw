import Box from "@mui/material/Box";
import AppBar from "@mui/material/AppBar";
import CssBaseline from "@mui/material/CssBaseline";
import Toolbar from "@mui/material/Toolbar";
import Typography from "@mui/material/Typography";
import IconButton from "@mui/material/IconButton";
import Tooltip from "@mui/material/Tooltip";
import Badge from "@mui/material/Badge";
import ChatIcon from "@mui/icons-material/Chat";
import CodeIcon from "@mui/icons-material/Code";
import LightModeIcon from "@mui/icons-material/LightMode";
import DarkModeIcon from "@mui/icons-material/DarkMode";
import BrushIcon from "@mui/icons-material/Brush";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import FrameProgressBar from "../components/ProgressBar";
import SideBar from "../components/SideBar";
import RoomManagementMenu from "../components/RoomManagementMenu";

import { useSocketManager } from "../hooks/useSocketManager";
import { useKeyboardShortcuts } from "../hooks/useKeyboardShortcuts";
import { useDragAndDrop } from "../hooks/useDragAndDrop";
import MyScene from "../components/Canvas";
import ChatWindow from "../components/ChatWindow";
import ConnectionDialog from "../components/ConnectionDialog";
import DropOverlay from "../components/DropOverlay";
import { useAppStore } from "../store";
import { useRestJoinManager } from "../hooks/useRestManager";
import React, { useState } from "react";
import { useColorScheme } from "@mui/material/styles";
import WindowManager from "../components/WindowManager";
import AddPlotButton from "../components/AddPlotButton";
import { useQueryClient } from "@tanstack/react-query";
import { LAYOUT_CONSTANTS } from "../constants/layout";

export default function MainPage() {
  useSocketManager();
  useKeyboardShortcuts();
  useRestJoinManager();
  const { isDragging, handleDragOver, handleDragEnter, handleDragLeave, handleDrop } = useDragAndDrop();

  const { chatOpen, setChatOpen, isDrawing, toggleDrawingMode, chatUnreadCount, serverVersion } = useAppStore();
  const [connectionDialogOpen, setConnectionDialogOpen] = useState(false);
  const { mode, setMode } = useColorScheme();
  const queryClient = useQueryClient();

  // File upload ref for button click
  const fileInputRef = React.useRef<HTMLInputElement>(null);

  const handleFileUploadClick = () => {
    fileInputRef.current?.click();
  };

  const handleFileInputChange = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const files = event.target.files;
    if (!files || files.length === 0) return;

    // Use the uploadFile API from drag/drop hook
    const file = files[0];
    const uploadEvent = {
      preventDefault: () => {},
      stopPropagation: () => {},
      dataTransfer: { files: [file] }
    } as any;

    await handleDrop(uploadEvent);

    // Reset input
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
    }
  };

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
          sx={{
            zIndex: (theme) => theme.zIndex.drawer + 1,
            height: LAYOUT_CONSTANTS.APPBAR_HEIGHT,
          }}
        >
          <Toolbar
            sx={{
              minHeight: `${LAYOUT_CONSTANTS.APPBAR_HEIGHT}px !important`,
              height: LAYOUT_CONSTANTS.APPBAR_HEIGHT,
              padding: '0 16px',
            }}
          >
            <Typography
              variant="h6"
              noWrap
              component="div"
              sx={{ flexGrow: 1 }}
            >
              {serverVersion ? `ZnDraw ${serverVersion}` : 'ZnDraw'}
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
                <Badge
                  badgeContent={chatUnreadCount}
                  color="error"
                  max={99}
                  invisible={chatUnreadCount === 0 || chatOpen}
                >
                  <ChatIcon />
                </Badge>
              </IconButton>
            </Tooltip>
            <AddPlotButton />
            <Tooltip title="Upload file">
              <IconButton
                color="inherit"
                aria-label="upload file"
                onClick={handleFileUploadClick}
              >
                <UploadFileIcon />
              </IconButton>
            </Tooltip>
            <RoomManagementMenu />
          </Toolbar>
        </AppBar>

        {/* Hidden file input for button upload */}
        <input
          type="file"
          ref={fileInputRef}
          style={{ display: 'none' }}
          onChange={handleFileInputChange}
          accept=".xyz,.extxyz,.pdb,.cif,.h5,.h5md,.hdf5,.gro,.mol,.sdf,.db,.json,.traj,.nc,.car,.xsf,.cube,.vasp,.poscar,.contcar,.xdatcar,.outcar,.xml,.pwi,.pwo,.out,.castep,.cell,.geom,.md,.gjf,.com,.log,.arc,.dmol"
        />

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
              onDragOver={handleDragOver}
              onDragEnter={handleDragEnter}
              onDragLeave={handleDragLeave}
              onDrop={handleDrop}
              sx={{
                flexGrow: 1,
                position: "relative",
                overflow: "hidden",
                display: "flex",
                flexDirection: "column",
              }}
            >
              <DropOverlay isDragging={isDragging} />
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
