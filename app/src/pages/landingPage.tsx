import Box from '@mui/material/Box';
import AppBar from '@mui/material/AppBar';
import CssBaseline from '@mui/material/CssBaseline';
import Toolbar from '@mui/material/Toolbar';
import Typography from '@mui/material/Typography';
import IconButton from '@mui/material/IconButton';
import Tooltip from '@mui/material/Tooltip';
import ChatIcon from '@mui/icons-material/Chat';
import CodeIcon from '@mui/icons-material/Code';
import LightModeIcon from '@mui/icons-material/LightMode';
import DarkModeIcon from '@mui/icons-material/DarkMode';
import FrameProgressBar from '../components/ProgressBar';
import SideBar from '../components/SideBar';

import { useSocketManager } from '../hooks/useSocketManager';
import { useKeyboardShortcuts } from '../hooks/useKeyboardShortcuts';
import MyScene from '../components/Canvas';
import ChatWindow from '../components/ChatWindow';
import ConnectionDialog from '../components/ConnectionDialog';
import { useAppStore } from '../store';
import { useRestJoinManager } from '../hooks/useRestManager';
import { useState } from 'react';
import { useColorScheme } from '@mui/material/styles';
import WindowManager from '../components/WindowManager';
import AddPlotButton from '../components/AddPlotButton';



export default function MainPage() {
  useSocketManager();
  useKeyboardShortcuts();
  useRestJoinManager();

  const { chatOpen, setChatOpen } = useAppStore();
  const [connectionDialogOpen, setConnectionDialogOpen] = useState(false);
  const { mode, setMode } = useColorScheme();

  const handleToggleColorMode = () => {
    setMode(mode === 'light' ? 'dark' : 'light');
  };

  return (
    <>
      <Box sx={{ display: 'flex' }}>
        <CssBaseline />
        <WindowManager />
        <AppBar position="fixed" sx={{ zIndex: (theme) => theme.zIndex.drawer + 1 }}>
          <Toolbar>
            <Typography variant="h6" noWrap component="div" sx={{ flexGrow: 1 }}>
              ZnDraw
            </Typography>
            <Tooltip title={mode === 'light' ? 'Switch to dark mode' : 'Switch to light mode'}>
              <IconButton
                color="inherit"
                aria-label="toggle color mode"
                onClick={handleToggleColorMode}
              >
                {mode === 'light' ? <DarkModeIcon /> : <LightModeIcon />}
              </IconButton>
            </Tooltip>
            <IconButton
              color="inherit"
              aria-label="show connection info"
              onClick={() => setConnectionDialogOpen(true)}
            >
              <CodeIcon />
            </IconButton>
            <IconButton
              color="inherit"
              aria-label="toggle chat"
              onClick={() => setChatOpen(!chatOpen)}
            >
              <ChatIcon />
            </IconButton>
            <AddPlotButton />
          </Toolbar>
        </AppBar>
        <SideBar />

        {/* Main Content Area */}
        <Box component="main" sx={{ flexGrow: 1, display: 'flex', flexDirection: 'column' }}>
          <Toolbar />
          <MyScene />
        </Box>

        {/* Chat Window */}
        <ChatWindow open={chatOpen} onClose={() => setChatOpen(false)} />
      </Box>
      <FrameProgressBar />

      <ConnectionDialog
        open={connectionDialogOpen}
        onClose={() => setConnectionDialogOpen(false)}
      />
    </>
  );
}