import React, { useState, useEffect } from "react";
import { useParams } from "react-router-dom";
import {
  Box,
  Container,
  Typography,
  TextField,
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  Card,
  CardContent,
  Grid,
  Chip,
  CircularProgress,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Paper,
  Tooltip,
} from "@mui/material";
import {
  Search as SearchIcon,
  CheckCircle,
  Error,
} from "@mui/icons-material";
import {
  useRoomExtensionsOverview,
  useGlobalExtensionsOverview,
} from "../hooks/useExtensionsOverview";
import { ExtensionDetailDrawer } from "../components/ExtensionDetailDrawer";

interface Props {
  mode: "room" | "global";
}

export const ExtensionsOverview: React.FC<Props> = ({ mode }) => {
  const { roomId } = useParams<{ roomId: string }>();

  const [search, setSearch] = useState("");
  const [category, setCategory] = useState<string>("all");
  const [selectedExtension, setSelectedExtension] = useState<{
    name: string;
    category: string;
  } | null>(null);

  // Debounced search
  const [debouncedSearch, setDebouncedSearch] = useState("");
  useEffect(() => {
    const timer = setTimeout(() => setDebouncedSearch(search), 300);
    return () => clearTimeout(timer);
  }, [search]);

  const filters = {
    category: category === "all" ? undefined : category,
    search: debouncedSearch || undefined,
  };

  const roomQuery = useRoomExtensionsOverview(roomId!, filters);
  const globalQuery = useGlobalExtensionsOverview(filters);

  const query = mode === "room" ? roomQuery : globalQuery;
  const data = query.data;

  if (query.isLoading) {
    return (
      <Box display="flex" justifyContent="center" alignItems="center" minHeight="400px">
        <CircularProgress />
      </Box>
    );
  }

  if (query.isError) {
    return (
      <Box p={4}>
        <Typography color="error">Failed to load extensions</Typography>
      </Box>
    );
  }

  return (
    <Container maxWidth="xl" sx={{ py: 4 }}>
      {/* Header */}
      <Box mb={4}>
        <Typography variant="h4" gutterBottom>
          {mode === "room" ? "Room Extensions" : "Global Extensions"}
        </Typography>
      </Box>

      {/* Filters */}
      <Box mb={3} display="flex" gap={2}>
        <TextField
          placeholder="Search extensions..."
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          InputProps={{
            startAdornment: <SearchIcon sx={{ mr: 1, color: "text.secondary" }} />,
          }}
          sx={{ flexGrow: 1 }}
        />
        <FormControl sx={{ minWidth: 200 }}>
          <InputLabel>Category</InputLabel>
          <Select
            value={category}
            onChange={(e) => setCategory(e.target.value)}
            label="Category"
          >
            <MenuItem value="all">All Categories</MenuItem>
            <MenuItem value="modifiers">Modifiers</MenuItem>
            <MenuItem value="selections">Selections</MenuItem>
            <MenuItem value="analysis">Analysis</MenuItem>
          </Select>
        </FormControl>
      </Box>

      {/* Summary Cards (room mode only) */}
      {mode === "room" && data && "summary" in data && (
        <Grid container spacing={2} mb={3}>
          <Grid item xs={12} sm={4}>
            <Card>
              <CardContent>
                <Typography color="text.secondary" gutterBottom>
                  Total Extensions
                </Typography>
                <Typography variant="h4">{data.summary.total_extensions}</Typography>
              </CardContent>
            </Card>
          </Grid>
          <Grid item xs={12} sm={4}>
            <Card>
              <CardContent>
                <Typography color="text.secondary" gutterBottom>
                  Active Workers
                </Typography>
                <Typography variant="h4">{data.summary.active_workers}</Typography>
              </CardContent>
            </Card>
          </Grid>
          <Grid item xs={12} sm={4}>
            <Card>
              <CardContent>
                <Typography color="text.secondary" gutterBottom>
                  Jobs (24h)
                </Typography>
                <Typography variant="h4">{data.summary.total_jobs_24h}</Typography>
              </CardContent>
            </Card>
          </Grid>
        </Grid>
      )}

      {/* Extensions Table */}
      <TableContainer component={Paper}>
        <Table>
          <TableHead>
            <TableRow>
              <TableCell>Name</TableCell>
              <TableCell>Category</TableCell>
              <TableCell>Provider</TableCell>
              {mode === "room" && <TableCell>Workers</TableCell>}
              {mode === "room" && <TableCell>Queue</TableCell>}
              {mode === "global" && <TableCell>Rooms</TableCell>}
              <TableCell>Total Jobs</TableCell>
              <TableCell>Success Rate</TableCell>
              {mode === "room" && <TableCell>Avg Wait</TableCell>}
              {mode === "room" && <TableCell>Avg Exec</TableCell>}
            </TableRow>
          </TableHead>
          <TableBody>
            {data?.extensions.map((ext: any) => (
              <TableRow
                key={`${ext.category}:${ext.name}`}
                hover
                onClick={() => setSelectedExtension({ name: ext.name, category: ext.category })}
                sx={{ cursor: "pointer" }}
              >
                <TableCell>{ext.name}</TableCell>
                <TableCell>
                  <Chip label={ext.category} size="small" />
                </TableCell>
                <TableCell>
                  <Chip
                    label={ext.provider === "celery" ? "Server" : "Client"}
                    size="small"
                    color={ext.provider === "celery" ? "primary" : "secondary"}
                  />
                </TableCell>
                {mode === "room" && (
                  <TableCell>
                    <Tooltip title="Idle workers / Progressing workers (actively executing jobs)">
                      <span>
                        {ext.workers.idle_count}/{ext.workers.progressing_count}
                      </span>
                    </Tooltip>
                  </TableCell>
                )}
                {mode === "room" && (
                  <TableCell>{ext.workers.queue_length}</TableCell>
                )}
                {mode === "global" && (
                  <TableCell>{ext.rooms.length}</TableCell>
                )}
                <TableCell>
                  {mode === "room" ? ext.analytics.total_jobs : ext.global_stats.total_jobs}
                </TableCell>
                <TableCell>
                  <Box display="flex" alignItems="center" gap={1}>
                    {(mode === "room"
                      ? ext.analytics.success_rate
                      : ext.global_stats.avg_success_rate
                    ).toFixed(1)}
                    %
                    {mode === "room" &&
                      (ext.analytics.success_rate > 90 ? (
                        <CheckCircle color="success" fontSize="small" />
                      ) : (
                        <Error color="error" fontSize="small" />
                      ))}
                  </Box>
                </TableCell>
                {mode === "room" && (
                  <TableCell>
                    {ext.analytics.avg_wait_time_ms > 0
                      ? `${(ext.analytics.avg_wait_time_ms / 1000).toFixed(2)}s`
                      : "N/A"}
                  </TableCell>
                )}
                {mode === "room" && (
                  <TableCell>
                    {ext.analytics.avg_execution_time_ms > 0
                      ? `${(ext.analytics.avg_execution_time_ms / 1000).toFixed(2)}s`
                      : "N/A"}
                  </TableCell>
                )}
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </TableContainer>

      {/* Detail Drawer */}
      {selectedExtension && mode === "room" && (
        <ExtensionDetailDrawer
          open={!!selectedExtension}
          onClose={() => setSelectedExtension(null)}
          roomId={roomId!}
          category={selectedExtension.category}
          extension={selectedExtension.name}
        />
      )}
    </Container>
  );
};
