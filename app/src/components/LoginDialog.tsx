import React, { useState } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  TextField,
  Button,
  Box,
  Typography,
  Alert,
} from '@mui/material';
import { login } from '../utils/auth';
import { useAppStore } from '../store';

interface LoginDialogProps {
  open: boolean;
  onClose: () => void;
}

export default function LoginDialog({ open, onClose }: LoginDialogProps) {
  const [username, setUsername] = useState('');
  const [password, setPassword] = useState('');
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  // Use individual selectors to prevent unnecessary re-renders
  const setUserName = useAppStore((state) => state.setUserName);
  const setUserRole = useAppStore((state) => state.setUserRole);
  const showSnackbar = useAppStore((state) => state.showSnackbar);

  const handleLoginAsGuest = async () => {
    setError(null);
    setLoading(true);
    try {
      // Auto-login without password (guest mode)
      const response = await login();
      setUserName(response.userName);
      setUserRole(response.role);

      // Socket will reconnect automatically when userName changes (via useSocketManager)

      showSnackbar(`Logged in as ${response.userName}`, 'success');
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Login failed');
    } finally {
      setLoading(false);
    }
  };

  const handleLogin = async () => {
    if (!username.trim()) {
      setError('Username is required');
      return;
    }
    if (!password) {
      setError('Password is required');
      return;
    }

    setError(null);
    setLoading(true);
    try {
      const response = await login(username, password);
      setUserName(response.userName);
      setUserRole(response.role);

      // Socket will reconnect automatically when userName changes (via useSocketManager)

      showSnackbar(
        response.role === 'admin' ? 'Logged in as admin' : 'Logged in successfully',
        'success'
      );
      onClose();

      // Clear form
      setUsername('');
      setPassword('');
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Login failed');
    } finally {
      setLoading(false);
    }
  };

  const handleKeyPress = (event: React.KeyboardEvent) => {
    if (event.key === 'Enter' && !loading) {
      if (username && password) {
        handleLogin();
      }
    }
  };

  const handleClose = () => {
    if (!loading) {
      setError(null);
      setUsername('');
      setPassword('');
      onClose();
    }
  };

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="xs" fullWidth>
      <DialogTitle>Login to ZnDraw</DialogTitle>
      <DialogContent>
        <Box sx={{ pt: 1, display: 'flex', flexDirection: 'column', gap: 2 }}>
          <Typography variant="body2" color="text.secondary">
            Login with your username and password, or continue as guest.
          </Typography>

          {error && (
            <Alert severity="error" onClose={() => setError(null)}>
              {error}
            </Alert>
          )}

          <TextField
            label="Username"
            value={username}
            onChange={(e) => setUsername(e.target.value)}
            onKeyPress={handleKeyPress}
            disabled={loading}
            fullWidth
            autoComplete="username"
          />

          <TextField
            label="Password"
            type="password"
            value={password}
            onChange={(e) => setPassword(e.target.value)}
            onKeyPress={handleKeyPress}
            disabled={loading}
            fullWidth
            autoComplete="current-password"
          />
        </Box>
      </DialogContent>
      <DialogActions sx={{ px: 3, pb: 2, flexDirection: 'column', gap: 1 }}>
        <Button
          onClick={handleLogin}
          disabled={loading || !username || !password}
          variant="contained"
          fullWidth
        >
          {loading ? 'Logging in...' : 'Login'}
        </Button>
        <Button
          onClick={handleLoginAsGuest}
          disabled={loading}
          variant="outlined"
          fullWidth
        >
          Continue as Guest
        </Button>
      </DialogActions>
    </Dialog>
  );
}
