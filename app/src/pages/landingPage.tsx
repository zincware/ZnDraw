import Box from "@mui/material/Box";
import AppBar from "@mui/material/AppBar";
import CssBaseline from "@mui/material/CssBaseline";
import Toolbar from "@mui/material/Toolbar";
import Typography from "@mui/material/Typography";
import IconButton from "@mui/material/IconButton";
import Tooltip from "@mui/material/Tooltip";
import Badge from "@mui/material/Badge";
import Menu from "@mui/material/Menu";
import MenuItem from "@mui/material/MenuItem";
import ListItemIcon from "@mui/material/ListItemIcon";
import ListItemText from "@mui/material/ListItemText";
import ChatIcon from "@mui/icons-material/Chat";
import CodeIcon from "@mui/icons-material/Code";
import LightModeIcon from "@mui/icons-material/LightMode";
import DarkModeIcon from "@mui/icons-material/DarkMode";
import BrushIcon from "@mui/icons-material/Brush";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import AccountCircleIcon from "@mui/icons-material/AccountCircle";
import LoginIcon from "@mui/icons-material/Login";
import LogoutIcon from "@mui/icons-material/Logout";
import AdminPanelSettingsIcon from "@mui/icons-material/AdminPanelSettings";
import ManageAccountsIcon from "@mui/icons-material/ManageAccounts";
import FrameProgressBar from "../components/ProgressBar";
import SideBar from "../components/SideBar";
import RoomManagementMenu from "../components/RoomManagementMenu";
import LoginDialog from "../components/LoginDialog";
import RegisterDialog from "../components/RegisterDialog";
import AdminPanel from "../components/AdminPanel";
import UserProfileDialog from "../components/UserProfileDialog";
import SiMGenTutorialDialog from "../components/SiMGenTutorialDialog";
import SiMGenButtons from "../components/SiMGenButtons";

import { useSocketManager } from "../hooks/useSocketManager";
import { useKeyboardShortcuts } from "../hooks/useKeyboardShortcuts";
import { useDragAndDrop } from "../hooks/useDragAndDrop";
import MyScene from "../components/Canvas";
import ChatWindow from "../components/ChatWindow";
import ConnectionDialog from "../components/ConnectionDialog";
import DropOverlay from "../components/DropOverlay";
import ProgressNotifications from "../components/ProgressNotifications";
import { useAppStore } from "../store";
import { useRestJoinManager } from "../hooks/useRestManager";
import React, { useState, useEffect } from "react";
import { useColorScheme } from "@mui/material/styles";
import { useLocation, useNavigate } from "react-router-dom";
import WindowManager from "../components/WindowManager";
import AddPlotButton from "../components/AddPlotButton";
import { useQueryClient, useMutation } from "@tanstack/react-query";
import { LAYOUT_CONSTANTS } from "../constants/layout";
import { getUsername, logout as authLogout, login as authLogin, getUserRole } from "../utils/auth";
import { loadFilesystemFile, loadGlobalFilesystemFile } from "../myapi/client";
import Link from "@mui/material/Link";

export default function MainPage() {
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const location = useLocation();
  const navigate = useNavigate();

  useSocketManager({ roomId: roomId || undefined });
  useKeyboardShortcuts();
  useRestJoinManager();
  const { isDragging, handleDragOver, handleDragEnter, handleDragLeave, handleDrop } = useDragAndDrop();

  const chatOpen = useAppStore((state) => state.chatOpen);
  const setChatOpen = useAppStore((state) => state.setChatOpen);
  const interactionMode = useAppStore((state) => state.mode);
  const enterDrawingMode = useAppStore((state) => state.enterDrawingMode);
  const exitDrawingMode = useAppStore((state) => state.exitDrawingMode);
  const chatUnreadCount = useAppStore((state) => state.chatUnreadCount);
  const serverVersion = useAppStore((state) => state.serverVersion);
  const globalSettings = useAppStore((state) => state.globalSettings);
  const setUserName = useAppStore((state) => state.setUserName);
  const setUserRole = useAppStore((state) => state.setUserRole);
  const showSnackbar = useAppStore((state) => state.showSnackbar);
  const [connectionDialogOpen, setConnectionDialogOpen] = useState(false);
  const [loginDialogOpen, setLoginDialogOpen] = useState(false);
  const [registerDialogOpen, setRegisterDialogOpen] = useState(false);
  const [adminPanelOpen, setAdminPanelOpen] = useState(false);
  const [userProfileDialogOpen, setUserProfileDialogOpen] = useState(false);
  const [profileAnchorEl, setProfileAnchorEl] = useState<null | HTMLElement>(null);
  const profileMenuOpen = Boolean(profileAnchorEl);
  const { mode: colorMode, setMode: setColorMode } = useColorScheme();
  const queryClient = useQueryClient();

  // Tutorial dialog state
  const [tutorialDialogOpen, setTutorialDialogOpen] = useState(false);

  // Filesystem load state
  const [pendingFilesystemLoad, setPendingFilesystemLoad] = useState<any>(null);

  // Mutation for loading files from filesystem
  const filesystemLoadMutation = useMutation({
    mutationFn: (data: { fsName: string; request: any; isPublic?: boolean }) => {
      // Route to correct endpoint based on isPublic flag
      return data.isPublic
        ? loadGlobalFilesystemFile(data.fsName, data.request)
        : loadFilesystemFile(roomId!, data.fsName, data.request);
    },
    onSuccess: (data) => {
      showSnackbar(`File loaded successfully: ${data.frameCount} frames`, "success");
    },
    onError: (error: any) => {
      showSnackbar(
        error?.response?.data?.error || "Failed to load file from filesystem",
        "error"
      );
    },
  });

  // Extract pending load from navigation state once
  useEffect(() => {
    const state = location.state as any;
    if (state?.pendingFilesystemLoad) {
      setPendingFilesystemLoad(state.pendingFilesystemLoad);
      // Clear navigation state immediately
      navigate(location.pathname, { replace: true, state: {} });
    }
  }, [location.pathname]); // Only when pathname changes (new navigation)

  // Trigger load when we have both roomId and pending load
  useEffect(() => {
    if (pendingFilesystemLoad && roomId) {
      const { fsName, request, isPublic } = pendingFilesystemLoad;

      // Clear pending load so this only runs once
      setPendingFilesystemLoad(null);

      // Trigger the mutation with isPublic flag
      filesystemLoadMutation.mutate({ fsName, request, isPublic });
    }
  }, [pendingFilesystemLoad, roomId]);

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

    handleDrop(uploadEvent);

    // Reset input
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
    }
  };

  const handleToggleColorMode = () => {
    setColorMode(colorMode === "light" ? "dark" : "light");
  };

  const handleProfileClick = (event: React.MouseEvent<HTMLElement>) => {
    setProfileAnchorEl(event.currentTarget);
  };

  const handleProfileClose = () => {
    setProfileAnchorEl(null);
  };

  const handleSwitchUser = () => {
    handleProfileClose();
    setLoginDialogOpen(true);
  };

  const handleLogout = async () => {
    handleProfileClose();

    // Auto-login as guest first
    try {
      authLogout();
      const response = await authLogin();

      // Update store - this will trigger useSocketManager's useEffect to reconnect
      setUserName(response.userName);
      setUserRole(response.role);

      // The socket will reconnect automatically when userName changes
      // The useSocketManager hook watches userName and will cleanup/reconnect

      showSnackbar('Logged out, reconnected as guest', 'info');
    } catch (err) {
      showSnackbar('Logout failed', 'error');
    }
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
            >
              {serverVersion ? `ZnDraw ${serverVersion}` : 'ZnDraw'}
              {globalSettings?.simgen?.enabled && (
                <>
                  {' + '}
                  <Link
                    href="https://github.com/RokasEl/simgen"
                    target="_blank"
                    rel="noopener noreferrer"
                    color="inherit"
                    underline="hover"
                    sx={{ cursor: 'pointer' }}
                  >
                    SiMGen
                  </Link>
                </>
              )}
            </Typography>
            <Box sx={{ flexGrow: 1 }} />
            <SiMGenButtons onTutorialClick={() => setTutorialDialogOpen(true)} />
            <Box sx={{ flexGrow: 1 }} />
            <Tooltip title={interactionMode === 'drawing' ? "Disable drawing mode" : "Enable drawing mode"}>
              <IconButton
                color="inherit"
                aria-label="toggle drawing mode"
                onClick={() => interactionMode === 'drawing' ? exitDrawingMode() : enterDrawingMode(queryClient)}
                sx={{
                  backgroundColor: interactionMode === 'drawing' ? 'rgba(255, 255, 255, 0.2)' : 'transparent',
                }}
              >
                <BrushIcon />
              </IconButton>
            </Tooltip>
            <Tooltip
              title={
                colorMode === "light"
                  ? "Switch to dark mode"
                  : "Switch to light mode"
              }
            >
              <IconButton
                color="inherit"
                aria-label="toggle color mode"
                onClick={handleToggleColorMode}
              >
                {colorMode === "light" ? <DarkModeIcon /> : <LightModeIcon />}
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
            <Tooltip title="User profile">
              <IconButton
                color="inherit"
                aria-label="user profile"
                onClick={handleProfileClick}
              >
                <AccountCircleIcon />
              </IconButton>
            </Tooltip>
          </Toolbar>
        </AppBar>

        {/* Profile Menu */}
        <Menu
          anchorEl={profileAnchorEl}
          open={profileMenuOpen}
          onClose={handleProfileClose}
          anchorOrigin={{
            vertical: 'bottom',
            horizontal: 'right',
          }}
          transformOrigin={{
            vertical: 'top',
            horizontal: 'right',
          }}
        >
          <MenuItem disabled>
            <ListItemText
              primary={getUsername() || 'Guest'}
              secondary={getUserRole() === 'admin' ? 'Admin' : getUserRole() === 'user' ? 'User' : 'Guest'}
            />
          </MenuItem>

          {/* Guest user menu items */}
          {getUserRole() === 'guest' && (
            <>
              <MenuItem onClick={() => { handleProfileClose(); setRegisterDialogOpen(true); }}>
                <ListItemIcon>
                  <ManageAccountsIcon fontSize="small" />
                </ListItemIcon>
                <ListItemText>Register</ListItemText>
              </MenuItem>
              <MenuItem onClick={handleSwitchUser}>
                <ListItemIcon>
                  <LoginIcon fontSize="small" />
                </ListItemIcon>
                <ListItemText>Login</ListItemText>
              </MenuItem>
            </>
          )}

          {/* Registered user menu items */}
          {getUserRole() !== 'guest' && (
            <MenuItem onClick={() => { handleProfileClose(); setUserProfileDialogOpen(true); }}>
              <ListItemIcon>
                <ManageAccountsIcon fontSize="small" />
              </ListItemIcon>
              <ListItemText>Manage Account</ListItemText>
            </MenuItem>
          )}

          {/* Admin-specific menu item */}
          {getUserRole() === 'admin' && (
            <MenuItem onClick={() => { handleProfileClose(); setAdminPanelOpen(true); }}>
              <ListItemIcon>
                <AdminPanelSettingsIcon fontSize="small" />
              </ListItemIcon>
              <ListItemText>Admin Panel</ListItemText>
            </MenuItem>
          )}

          {/* Logout - available to all */}
          <MenuItem onClick={handleLogout}>
            <ListItemIcon>
              <LogoutIcon fontSize="small" />
            </ListItemIcon>
            <ListItemText>Logout</ListItemText>
          </MenuItem>
        </Menu>

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

      <LoginDialog
        open={loginDialogOpen}
        onClose={() => setLoginDialogOpen(false)}
      />

      <RegisterDialog
        open={registerDialogOpen}
        onClose={() => setRegisterDialogOpen(false)}
      />

      <AdminPanel
        open={adminPanelOpen}
        onClose={() => setAdminPanelOpen(false)}
      />

      <UserProfileDialog
        open={userProfileDialogOpen}
        onClose={() => setUserProfileDialogOpen(false)}
      />

      <SiMGenTutorialDialog
        open={tutorialDialogOpen}
        onClose={() => setTutorialDialogOpen(false)}
        url="https://slides.com/rokasel/zndrawtutorial-9cc179/fullscreen?style=light"
      />

      {/* Progress notifications */}
      <ProgressNotifications />
    </>
  );
}
