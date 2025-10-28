import React from "react";
import {
  Drawer,
  Box,
  Typography,
  IconButton,
  Divider,
  Chip,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  CircularProgress,
} from "@mui/material";
import { Close as CloseIcon } from "@mui/icons-material";
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
} from "recharts";
import { useExtensionDetailedAnalytics } from "../hooks/useExtensionsOverview";

interface Props {
  open: boolean;
  onClose: () => void;
  roomId: string;
  category: string;
  extension: string;
}

export const ExtensionDetailDrawer: React.FC<Props> = ({
  open,
  onClose,
  roomId,
  category,
  extension,
}) => {
  const { data, isLoading } = useExtensionDetailedAnalytics(
    roomId,
    category,
    extension
  );

  return (
    <Drawer anchor="right" open={open} onClose={onClose}>
      <Box sx={{ width: 600, p: 3 }}>
        {/* Header */}
        <Box display="flex" justifyContent="space-between" alignItems="center" mb={2}>
          <Box>
            <Typography variant="h5">{extension}</Typography>
            <Chip label={category} size="small" sx={{ mt: 1 }} />
          </Box>
          <IconButton onClick={onClose}>
            <CloseIcon />
          </IconButton>
        </Box>

        <Divider sx={{ mb: 3 }} />

        {isLoading ? (
          <Box display="flex" justifyContent="center" py={4}>
            <CircularProgress />
          </Box>
        ) : (
          <>
            {/* Summary Stats */}
            <Box mb={3}>
              <Typography variant="h6" gutterBottom>
                Statistics
              </Typography>
              <Box display="grid" gridTemplateColumns="1fr 1fr" gap={2}>
                <Box>
                  <Typography color="text.secondary" variant="body2">
                    Total Jobs
                  </Typography>
                  <Typography variant="h6">
                    {data?.total_stats.total_jobs || 0}
                  </Typography>
                </Box>
                <Box>
                  <Typography color="text.secondary" variant="body2">
                    Success Rate
                  </Typography>
                  <Typography variant="h6">
                    {data?.total_stats.overall_success_rate.toFixed(1)}%
                  </Typography>
                </Box>
                <Box>
                  <Typography color="text.secondary" variant="body2">
                    Avg Wait Time
                  </Typography>
                  <Typography variant="h6">
                    {data?.total_stats.overall_avg_wait_ms
                      ? `${(data.total_stats.overall_avg_wait_ms / 1000).toFixed(2)}s`
                      : "N/A"}
                  </Typography>
                </Box>
                <Box>
                  <Typography color="text.secondary" variant="body2">
                    Avg Execution Time
                  </Typography>
                  <Typography variant="h6">
                    {data?.total_stats.overall_avg_execution_ms
                      ? `${(data.total_stats.overall_avg_execution_ms / 1000).toFixed(2)}s`
                      : "N/A"}
                  </Typography>
                </Box>
              </Box>
            </Box>

            {/* Usage Chart */}
            <Box mb={3}>
              <Typography variant="h6" gutterBottom>
                Usage Trend
              </Typography>
              <ResponsiveContainer width="100%" height={200}>
                <LineChart data={data?.daily_stats || []}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis dataKey="date" />
                  <YAxis />
                  <Tooltip />
                  <Line
                    type="monotone"
                    dataKey="job_count"
                    stroke="#8884d8"
                    name="Jobs"
                  />
                </LineChart>
              </ResponsiveContainer>
            </Box>

            {/* Error Breakdown */}
            {data?.error_breakdown && data.error_breakdown.length > 0 && (
              <Box mb={3}>
                <Typography variant="h6" gutterBottom>
                  Recent Errors
                </Typography>
                <TableContainer>
                  <Table size="small">
                    <TableHead>
                      <TableRow>
                        <TableCell>Error</TableCell>
                        <TableCell align="right">Count</TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {data.error_breakdown.slice(0, 5).map((err, idx) => (
                        <TableRow key={idx}>
                          <TableCell>{err.error}</TableCell>
                          <TableCell align="right">{err.count}</TableCell>
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>
                </TableContainer>
              </Box>
            )}
          </>
        )}
      </Box>
    </Drawer>
  );
};
