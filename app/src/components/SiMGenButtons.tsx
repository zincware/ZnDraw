import { useEffect, useState } from 'react';
import Box from '@mui/material/Box';
import Tooltip from '@mui/material/Tooltip';
import Button from '@mui/material/Button';
import SchoolIcon from '@mui/icons-material/School';
import LinkIcon from '@mui/icons-material/Link';
import ScienceIcon from '@mui/icons-material/Science';
import ArticleIcon from '@mui/icons-material/Article';
import { submitExtension } from '../myapi/client';
import { useAppStore } from '../store';

interface SiMGenButtonsProps {
  onTutorialClick: () => void;
}

export default function SiMGenButtons({ onTutorialClick }: SiMGenButtonsProps) {
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const globalSettings = useAppStore((state) => state.globalSettings);
  const showSnackbar = useAppStore((state) => state.showSnackbar);
  const [showAnimation, setShowAnimation] = useState(true);

  // Disable animation after it plays once
  useEffect(() => {
    if (showAnimation) {
      const timer = setTimeout(() => {
        setShowAnimation(false);
      }, 3000); // Animation duration (3 pulses of 1s each)
      return () => clearTimeout(timer);
    }
  }, [showAnimation]);

  // Handler for submitting SiMGen extensions
  const handleExtensionSubmit = async (extensionName: string) => {
    if (!roomId) {
      showSnackbar('Please select a room first', 'warning');
      return;
    }

    try {
      // SiMGen extensions are registered publicly
      await submitExtension({
        roomId,
        category: 'modifiers',
        extension: extensionName,
        data: {},
        isPublic: true,
      });
      showSnackbar(`${extensionName} submitted successfully`, 'success');
    } catch (error) {
      console.error(`Failed to submit ${extensionName}:`, error);
      showSnackbar(`Failed to submit ${extensionName}`, 'error');
    }
  };

  // Don't render anything if simgen is not enabled
  if (!globalSettings?.simgen?.enabled) {
    return null;
  }

  return (
    <Box
      component="fieldset"
      sx={{
        display: 'flex',
        gap: 0.5,
        px: 2,
        py: 0,
        m: 0,
        borderRadius: 2,
        bgcolor: (theme) => theme.palette.mode === 'light'
          ? 'rgba(255, 191, 0, 0)'
          : 'rgba(255, 193, 7, 0.0)',
        border: (theme) => `1px solid ${
          theme.palette.mode === 'light'
            ? 'rgba(255, 193, 7, 0.3)'
            : 'rgba(255, 193, 7, 0.2)'
        }`,
        animation: showAnimation ? 'simgenPulse 1s ease-in-out 3' : 'none',
        '@keyframes simgenPulse': {
          '0%, 100%': {
            boxShadow: '0 0 0 0 rgba(255, 193, 7, 0)',
          },
          '50%': {
            boxShadow: '0 0 20px 4px rgba(255, 193, 7, 0.4)',
          },
        },
      }}
    >
      <Box
        component="legend"
        sx={{
          px: 1,
          m: 0,
          fontSize: '0.75rem',
          lineHeight: 1.2,
        }}
      >
        SiMGen
      </Box>
      <Tooltip title="SiMGen Tutorial">
        <Button
          color="inherit"
          startIcon={<SchoolIcon />}
          onClick={onTutorialClick}
          size="small"
        >
          Tutorial
        </Button>
      </Tooltip>
      <Tooltip title="Connect selected atoms">
        <Button
          color="inherit"
          startIcon={<LinkIcon />}
          onClick={() => handleExtensionSubmit('Connect')}
          size="small"
        >
          Connect
        </Button>
      </Tooltip>
      <Tooltip title="Run SiMGen Demo">
        <Button
          color="inherit"
          startIcon={<ScienceIcon />}
          onClick={() => handleExtensionSubmit('SiMGenDemo')}
          size="small"
        >
          SiMGen
        </Button>
      </Tooltip>
      <Tooltip title="Replace scene with empty canvas">
        <Button
          color="inherit"
          startIcon={<ArticleIcon />}
          onClick={() => handleExtensionSubmit('NewCanvas')}
          size="small"
        >
          New Canvas
        </Button>
      </Tooltip>
    </Box>
  );
}
