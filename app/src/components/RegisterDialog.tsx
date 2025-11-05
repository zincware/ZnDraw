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
import { registerUser, getUsername } from '../utils/auth';
import { useAppStore } from '../store';

interface RegisterDialogProps {
  open: boolean;
  onClose: () => void;
}

export default function RegisterDialog({ open, onClose }: RegisterDialogProps) {
  const [username, setUsername] = useState('');
  const [password, setPassword] = useState('');
  const [passwordConfirm, setPasswordConfirm] = useState('');
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  // Use individual selectors to prevent unnecessary re-renders
  const setUserName = useAppStore((state) => state.setUserName);
  const setUserRole = useAppStore((state) => state.setUserRole);
  const showSnackbar = useAppStore((state) => state.showSnackbar);

  const handleRegister = async () => {
    setError(null);

    if (!username || !username.trim()) {
      setError('Username is required');
      return;
    }

    if (!password) {
      setError('Password is required');
      return;
    }

    if (password !== passwordConfirm) {
      setError('Passwords do not match');
      return;
    }

    setLoading(true);
    try {
      const response = await registerUser(username, password);

      // Update store with new username and role
      setUserName(response.userName);
      setUserRole(response.role);

      // Socket will reconnect automatically when userName changes (via useSocketManager)

      showSnackbar(`Registered as ${response.userName}`, 'success');
      onClose();

      // Clear form
      setUsername('');
      setPassword('');
      setPasswordConfirm('');
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Registration failed');
    } finally {
      setLoading(false);
    }
  };

  const handleKeyPress = (event: React.KeyboardEvent) => {
    if (event.key === 'Enter' && !loading) {
      if (username && password && passwordConfirm) {
        handleRegister();
      }
    }
  };

  const handleClose = () => {
    if (!loading) {
      setError(null);
      setUsername('');
      setPassword('');
      setPasswordConfirm('');
      onClose();
    }
  };

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="xs" fullWidth>
      <DialogTitle>Register Account</DialogTitle>
      <DialogContent>
        <Box sx={{ pt: 1, display: 'flex', flexDirection: 'column', gap: 2 }}>
          <Typography variant="body2" color="text.secondary">
            Current temporary username: <strong>{getUsername()}</strong>
          </Typography>
          <Typography variant="body2" color="text.secondary">
            Choose a permanent username and password to register your account.
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
            helperText="Cannot be changed after registration"
          />

          <TextField
            label="Password"
            type="password"
            value={password}
            onChange={(e) => setPassword(e.target.value)}
            onKeyPress={handleKeyPress}
            disabled={loading}
            fullWidth
            autoComplete="new-password"
          />

          <TextField
            label="Confirm Password"
            type="password"
            value={passwordConfirm}
            onChange={(e) => setPasswordConfirm(e.target.value)}
            onKeyPress={handleKeyPress}
            disabled={loading}
            fullWidth
            autoComplete="new-password"
          />
        </Box>
      </DialogContent>
      <DialogActions sx={{ px: 3, pb: 2 }}>
        <Button onClick={handleClose} disabled={loading}>
          Cancel
        </Button>
        <Button
          onClick={handleRegister}
          disabled={loading || !username || !password || !passwordConfirm}
          variant="contained"
        >
          {loading ? 'Registering...' : 'Register'}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
