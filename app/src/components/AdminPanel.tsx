import React, { useState, useEffect } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Paper,
  IconButton,
  Chip,
  TextField,
  Box,
  Typography,
  Alert,
  Tooltip,
} from '@mui/material';
import DeleteIcon from '@mui/icons-material/Delete';
import SecurityIcon from '@mui/icons-material/Security';
import PersonIcon from '@mui/icons-material/Person';
import LockResetIcon from '@mui/icons-material/LockReset';
import { getToken } from '../utils/auth';

interface User {
  userName: string;
  role: 'guest' | 'user' | 'admin';
  createdAt?: string;
  lastLogin?: string;
}

interface AdminPanelProps {
  open: boolean;
  onClose: () => void;
}

export default function AdminPanel({ open, onClose }: AdminPanelProps) {
  const [users, setUsers] = useState<User[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [resetPasswordDialog, setResetPasswordDialog] = useState<{
    open: boolean;
    userName: string;
  } | null>(null);
  const [newPassword, setNewPassword] = useState('');

  const fetchUsers = async () => {
    setLoading(true);
    setError(null);
    try {
      const token = getToken();
      const response = await fetch('/api/admin/users', {
        headers: {
          Authorization: `Bearer ${token}`,
        },
      });

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({ error: response.statusText }));
        throw new Error(errorData.error || 'Failed to fetch users');
      }

      const data = await response.json();
      setUsers(data.users);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load users');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    if (open) {
      fetchUsers();
    }
  }, [open]);

  const handlePromote = async (userName: string) => {
    try {
      const token = getToken();
      const response = await fetch(`/api/admin/users/${userName}/promote`, {
        method: 'POST',
        headers: {
          Authorization: `Bearer ${token}`,
        },
      });

      if (!response.ok) {
        throw new Error('Failed to promote user');
      }

      await fetchUsers();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to promote user');
    }
  };

  const handleDemote = async (userName: string) => {
    try {
      const token = getToken();
      const response = await fetch(`/api/admin/users/${userName}/demote`, {
        method: 'POST',
        headers: {
          Authorization: `Bearer ${token}`,
        },
      });

      if (!response.ok) {
        throw new Error('Failed to demote user');
      }

      await fetchUsers();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to demote user');
    }
  };

  const handleDelete = async (userName: string) => {
    if (!confirm(`Are you sure you want to delete user "${userName}"? This cannot be undone.`)) {
      return;
    }

    try {
      const token = getToken();
      const response = await fetch(`/api/admin/users/${userName}`, {
        method: 'DELETE',
        headers: {
          Authorization: `Bearer ${token}`,
        },
      });

      if (!response.ok) {
        throw new Error('Failed to delete user');
      }

      await fetchUsers();
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to delete user');
    }
  };

  const handleResetPasswordOpen = (userName: string) => {
    setResetPasswordDialog({ open: true, userName });
    setNewPassword('');
  };

  const handleResetPasswordClose = () => {
    setResetPasswordDialog(null);
    setNewPassword('');
  };

  const handleResetPasswordSubmit = async () => {
    if (!resetPasswordDialog || !newPassword) {
      return;
    }

    try {
      const token = getToken();
      const response = await fetch(`/api/admin/users/${resetPasswordDialog.userName}/reset-password`, {
        method: 'POST',
        headers: {
          Authorization: `Bearer ${token}`,
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ newPassword }),
      });

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({ error: response.statusText }));
        throw new Error(errorData.error || 'Failed to reset password');
      }

      handleResetPasswordClose();
      setError(null);
      // Show success message
      alert(`Password reset successfully for ${resetPasswordDialog.userName}`);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to reset password');
    }
  };

  const getRoleChip = (role: string) => {
    switch (role) {
      case 'admin':
        return <Chip label="Admin" color="error" size="small" />;
      case 'user':
        return <Chip label="User" color="primary" size="small" />;
      case 'guest':
        return <Chip label="Guest" color="default" size="small" />;
      default:
        return <Chip label={role} size="small" />;
    }
  };

  const formatDate = (dateString?: string) => {
    if (!dateString) return 'N/A';
    try {
      return new Date(dateString).toLocaleString();
    } catch {
      return 'Invalid date';
    }
  };

  return (
    <>
      <Dialog open={open} onClose={onClose} maxWidth="lg" fullWidth>
        <DialogTitle>User Management</DialogTitle>
        <DialogContent>
          <Box sx={{ mt: 1 }}>
            {error && (
              <Alert severity="error" onClose={() => setError(null)} sx={{ mb: 2 }}>
                {error}
              </Alert>
            )}

            {loading ? (
              <Typography>Loading users...</Typography>
            ) : (
              <TableContainer component={Paper}>
                <Table>
                  <TableHead>
                    <TableRow>
                      <TableCell>Username</TableCell>
                      <TableCell>Role</TableCell>
                      <TableCell>Created</TableCell>
                      <TableCell>Last Login</TableCell>
                      <TableCell align="right">Actions</TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {users.map((user) => (
                      <TableRow key={user.userName}>
                        <TableCell>{user.userName}</TableCell>
                        <TableCell>{getRoleChip(user.role)}</TableCell>
                        <TableCell>
                          <Typography variant="caption">{formatDate(user.createdAt)}</Typography>
                        </TableCell>
                        <TableCell>
                          <Typography variant="caption">{formatDate(user.lastLogin)}</Typography>
                        </TableCell>
                        <TableCell align="right">
                          <Box sx={{ display: 'flex', gap: 0.5, justifyContent: 'flex-end' }}>
                            {user.role === 'admin' ? (
                              <Tooltip title="Demote to user">
                                <IconButton
                                  size="small"
                                  onClick={() => handleDemote(user.userName)}
                                  color="warning"
                                >
                                  <PersonIcon />
                                </IconButton>
                              </Tooltip>
                            ) : (
                              <Tooltip title="Promote to admin">
                                <IconButton
                                  size="small"
                                  onClick={() => handlePromote(user.userName)}
                                  color="primary"
                                >
                                  <SecurityIcon />
                                </IconButton>
                              </Tooltip>
                            )}

                            {user.role !== 'guest' && (
                              <Tooltip title="Reset password">
                                <IconButton
                                  size="small"
                                  onClick={() => handleResetPasswordOpen(user.userName)}
                                  color="info"
                                >
                                  <LockResetIcon />
                                </IconButton>
                              </Tooltip>
                            )}

                            <Tooltip title="Delete user">
                              <IconButton
                                size="small"
                                onClick={() => handleDelete(user.userName)}
                                color="error"
                              >
                                <DeleteIcon />
                              </IconButton>
                            </Tooltip>
                          </Box>
                        </TableCell>
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              </TableContainer>
            )}
          </Box>
        </DialogContent>
        <DialogActions>
          <Button onClick={fetchUsers} disabled={loading}>
            Refresh
          </Button>
          <Button onClick={onClose}>Close</Button>
        </DialogActions>
      </Dialog>

      {/* Reset Password Dialog */}
      <Dialog
        open={resetPasswordDialog?.open || false}
        onClose={handleResetPasswordClose}
        maxWidth="xs"
        fullWidth
      >
        <DialogTitle>Reset Password</DialogTitle>
        <DialogContent>
          <Box sx={{ pt: 1 }}>
            <Typography variant="body2" color="text.secondary" gutterBottom>
              Reset password for user: <strong>{resetPasswordDialog?.userName}</strong>
            </Typography>
            <TextField
              autoFocus
              margin="dense"
              label="New Password"
              type="password"
              fullWidth
              value={newPassword}
              onChange={(e) => setNewPassword(e.target.value)}
              onKeyPress={(e) => {
                if (e.key === 'Enter' && newPassword) {
                  handleResetPasswordSubmit();
                }
              }}
            />
          </Box>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleResetPasswordClose}>Cancel</Button>
          <Button
            onClick={handleResetPasswordSubmit}
            variant="contained"
            disabled={!newPassword}
          >
            Reset Password
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}
