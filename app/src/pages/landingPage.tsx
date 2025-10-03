import Box from '@mui/material/Box';
import AppBar from '@mui/material/AppBar';
import CssBaseline from '@mui/material/CssBaseline';
import Toolbar from '@mui/material/Toolbar';
import Typography from '@mui/material/Typography';
import IconButton from '@mui/material/IconButton';
import ChatIcon from '@mui/icons-material/Chat';
import CodeIcon from '@mui/icons-material/Code';
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



export default function MainPage() {
  useSocketManager();
  useKeyboardShortcuts();
  useRestJoinManager();

  const { chatOpen, setChatOpen } = useAppStore();
  const [connectionDialogOpen, setConnectionDialogOpen] = useState(false);

  return (
    <>
      <Box sx={{ display: 'flex' }}>

        <CssBaseline />
        <AppBar position="fixed" sx={{ zIndex: (theme) => theme.zIndex.drawer + 1 }}>
          <Toolbar>
            <Typography variant="h6" noWrap component="div" sx={{ flexGrow: 1 }}>
              ZnDraw
            </Typography>
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