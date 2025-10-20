/**
 * PropertySelector component for selecting properties to inspect.
 * Displays categorized lists of per-particle and global properties with search.
 */

import { useState, useMemo } from "react";
import {
  Box,
  List,
  ListItem,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  Checkbox,
  TextField,
  Chip,
  Typography,
  Collapse,
} from "@mui/material";
import { ExpandMore, ExpandLess, Search as SearchIcon, Clear as ClearIcon } from "@mui/icons-material";
import type { PropertyInfo } from "../../../types/property-inspector";

interface PropertySelectorProps {
  perParticleProps: PropertyInfo[];
  globalProps: PropertyInfo[];
  selectedKeys: string[];
  onToggle: (key: string) => void;
  onClearAll: () => void;
}

export default function PropertySelector({
  perParticleProps,
  globalProps,
  selectedKeys,
  onToggle,
  onClearAll,
}: PropertySelectorProps) {
  const [searchQuery, setSearchQuery] = useState("");
  const [perParticleExpanded, setPerParticleExpanded] = useState(true);
  const [globalExpanded, setGlobalExpanded] = useState(true);

  const filteredPerParticle = useMemo(() => {
    if (!searchQuery) return perParticleProps;
    const query = searchQuery.toLowerCase();
    return perParticleProps.filter((p) => p.key.toLowerCase().includes(query));
  }, [perParticleProps, searchQuery]);

  const filteredGlobal = useMemo(() => {
    if (!searchQuery) return globalProps;
    const query = searchQuery.toLowerCase();
    return globalProps.filter((p) => p.key.toLowerCase().includes(query));
  }, [globalProps, searchQuery]);

  return (
    <Box sx={{ p: 1 }}>
      <TextField
        fullWidth
        size="small"
        placeholder="Search properties..."
        value={searchQuery}
        onChange={(e) => setSearchQuery(e.target.value)}
        slotProps={{
          input: {
            startAdornment: <SearchIcon sx={{ mr: 1, color: "text.secondary" }} />,
          },
        }}
        sx={{ mb: 1 }}
      />

      <Box sx={{ mb: 1, display: "flex", gap: 1, flexWrap: "wrap" }}>
        <Chip label={`${selectedKeys.length} selected`} size="small" color="primary" />
        <Chip label={`${perParticleProps.length} per-particle`} size="small" />
        <Chip label={`${globalProps.length} global`} size="small" />
        {selectedKeys.length > 0 && (
          <Chip
            label="Clear all"
            size="small"
            onDelete={onClearAll}
            deleteIcon={<ClearIcon />}
          />
        )}
      </Box>

      <List dense>
        <ListItem>
          <ListItemButton onClick={() => setPerParticleExpanded(!perParticleExpanded)}>
            <ListItemText
              primary={<Typography variant="subtitle2">Per-Particle Properties</Typography>}
            />
            {perParticleExpanded ? <ExpandLess /> : <ExpandMore />}
          </ListItemButton>
        </ListItem>
        <Collapse in={perParticleExpanded} timeout="auto" unmountOnExit>
          {filteredPerParticle.map((prop) => (
            <ListItem key={prop.key} dense>
              <ListItemButton onClick={() => onToggle(prop.key)}>
                <ListItemIcon>
                  <Checkbox
                    edge="start"
                    checked={selectedKeys.includes(prop.key)}
                    tabIndex={-1}
                    disableRipple
                  />
                </ListItemIcon>
                <ListItemText
                  primary={prop.key}
                  secondary={`${prop.metadata.dtype} [${prop.metadata.shape?.join(" × ") ?? "scalar"}]`}
                  slotProps={{
                    primary: { sx: { fontFamily: "monospace", fontSize: "0.85rem" } },
                  }}
                />
              </ListItemButton>
            </ListItem>
          ))}
        </Collapse>

        <ListItem>
          <ListItemButton onClick={() => setGlobalExpanded(!globalExpanded)}>
            <ListItemText
              primary={<Typography variant="subtitle2">Global Properties</Typography>}
            />
            {globalExpanded ? <ExpandLess /> : <ExpandMore />}
          </ListItemButton>
        </ListItem>
        <Collapse in={globalExpanded} timeout="auto" unmountOnExit>
          {filteredGlobal.map((prop) => (
            <ListItem key={prop.key} dense>
              <ListItemButton onClick={() => onToggle(prop.key)}>
                <ListItemIcon>
                  <Checkbox
                    edge="start"
                    checked={selectedKeys.includes(prop.key)}
                    tabIndex={-1}
                    disableRipple
                  />
                </ListItemIcon>
                <ListItemText
                  primary={prop.key}
                  secondary={`${prop.metadata.dtype} [${prop.metadata.shape?.join(" × ") ?? "scalar"}]`}
                  slotProps={{
                    primary: { sx: { fontFamily: "monospace", fontSize: "0.85rem" } },
                  }}
                />
              </ListItemButton>
            </ListItem>
          ))}
        </Collapse>
      </List>
    </Box>
  );
}
