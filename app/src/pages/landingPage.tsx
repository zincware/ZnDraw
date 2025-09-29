import Box from '@mui/material/Box';
import AppBar from '@mui/material/AppBar';
import CssBaseline from '@mui/material/CssBaseline';
import Toolbar from '@mui/material/Toolbar';
import Typography from '@mui/material/Typography';
import FrameProgressBar from '../components/ProgressBar';
import SideBar from '../components/SideBar';

import { useSocketManager } from '../hooks/useSocketManager';
import MyScene from '../components/Canvas';



export default function MainPage() {
  useSocketManager();

  return (
    <>
      <Box sx={{ display: 'flex' }}>

        <CssBaseline />
        <AppBar position="fixed" sx={{ zIndex: (theme) => theme.zIndex.drawer + 1 }}>
          <Toolbar>
            <Typography variant="h6" noWrap component="div">
              ZnDraw
            </Typography>
          </Toolbar>
        </AppBar>
        <SideBar />

        {/* Main Content Area */}
        <Box component="main" sx={{ flexGrow: 1, display: 'flex', flexDirection: 'column' }}>
          <Toolbar />
          <MyScene />
        </Box>
      </Box>
      <FrameProgressBar />
    </>
  );
}