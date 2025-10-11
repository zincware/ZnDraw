# Property Inspector Panel - Implementation Plan

## Overview

A custom JSONForms renderer integrated into the **Settings panel** for inspecting per-particle and global frame properties. Uses MUI-X DataGrid for efficient display with property toggles. Selected properties are fetched iteratively using React Query for optimal performance.

## Key Architecture Decision

**Property Inspector is NOT a separate sidebar panel.** Instead:
- ✅ `PropertyInspectorSettings` is part of `RoomConfig` Pydantic model
- ✅ Appears in the **Settings sidebar panel** via JSONForms
- ✅ Uses `x-custom-type="property-inspector"` custom renderer
- ✅ Follows existing patterns: `x-custom-type="dynamic-enum"`, `format="color"`, `format="range"`

## Problem Analysis

### Current Limitations
1. **Hover InfoBox** is mouse-bound - users can't interact with it
2. **No discoverability** - users don't know what properties are available
3. **Performance concerns** - fetching all properties would be expensive
4. **Type detection needed** - distinguish per-particle vs. global properties

### Data Structure (from API)

```json
{
  "frameId": 1,
  "keys": ["arrays.positions", "calc.forces", "cell", "pbc", ...],
  "metadata": {
    "arrays.positions": {
      "dtype": "float64",
      "shape": [3000, 3],  // First dim = particle count → per-particle
      "type": "array"
    },
    "calc.forces": {
      "dtype": "float64",
      "shape": [3000, 3],  // First dim = particle count → per-particle
      "type": "array"
    },
    "cell": {
      "dtype": "float64",
      "shape": [3, 3],     // First dim ≠ particle count → global
      "type": "array"
    },
    "pbc": {
      "dtype": "bool",
      "shape": [3],        // First dim ≠ particle count → global
      "type": "array"
    }
  }
}
```

**Detection Logic:**
```typescript
const isPerParticleProperty = (
  metadata: PropertyMetadata,
  particleCount: number
): boolean => {
  // Check if first dimension of shape equals particle count
  return metadata.shape[0] === particleCount;
};
```

---

## Architecture Design

### 1. Pydantic Settings Integration (Backend)

**Add PropertyInspector model to settings:**

```python
# src/zndraw/settings.py

class PropertyInspector(SettingsBase):
    """Property Inspector settings for per-particle and global property display."""
    
    enabled_properties: list[str] = Field(
        default_factory=list,
        description="Selected property keys to display in the inspector table"
    )
    
    show_per_particle: bool = Field(
        default=True,
        description="Show per-particle properties (arrays.*, calc.*)"
    )
    
    @classmethod
    def model_json_schema(cls, *args, **kwargs) -> dict[str, t.Any]:
        schema = super().model_json_schema(*args, **kwargs)
        # Mark enabled_properties for custom renderer
        schema["properties"]["enabled_properties"]["x-custom-type"] = "property-inspector"
        return schema


class RoomConfig(SettingsBase):
    """ZnDraw room configuration combining all settings sections."""
    
    camera: Camera = Camera()
    studio_lighting: StudioLighting = StudioLighting()
    property_inspector: PropertyInspector = PropertyInspector()  # NEW
```

**Key Points:**
- `enabled_properties` gets `x-custom-type="property-inspector"` in its JSON schema
- This triggers the custom PropertyInspectorRenderer in the frontend
- Other fields (`show_per_particle`, `show_global`) use standard boolean renderers
- Rendered automatically in Settings panel via JSONForms

---

### 2. Frontend Custom Renderer (JSONForms Pattern)

**Create custom JSONForms renderer:**

```
app/src/components/jsonforms-renderers/
├── PropertyInspectorRenderer.tsx     # Custom renderer for x-custom-type="property-inspector"
│   ├── PropertyInspectorTester       # Tester function (rankWith + schemaMatches)
│   └── PropertyInspectorRenderer     # React component (wrapped with withJsonFormsControlProps)
```

**Structure follows existing patterns:**
- Similar to `DynamicEnumRenderer.tsx` (ControlProps-based)
- Similar to `CustomColorPicker.tsx` (format="color")
- Similar to `CustomRangeSlider.tsx` (format="range")

---

## Implementation Details

### Phase 1: Backend - Pydantic Settings Model

#### 1.1 Create PropertyInspector Settings Class

```python
# src/zndraw/settings.py

class PropertyInspector(SettingsBase):
    """Property Inspector settings for per-particle and global property display."""
    
    enabled_properties: list[str] = Field(
        default_factory=list,
        description="Selected property keys to display in the inspector table"
    )
    
    show_per_particle: bool = Field(
        default=True,
        description="Show per-particle properties (shape[0] == particle_count)"
    )
    
    
    @classmethod
    def model_json_schema(cls, *args, **kwargs) -> dict[str, t.Any]:
        """Inject custom type for PropertyInspectorRenderer."""
        schema = super().model_json_schema(*args, **kwargs)
        # Mark enabled_properties field for custom renderer
        schema["properties"]["enabled_properties"]["x-custom-type"] = "property-inspector"
        return schema


# Add to RoomConfig
class RoomConfig(SettingsBase):
    """ZnDraw room configuration combining all settings sections."""
    
    camera: Camera = Camera()
    studio_lighting: StudioLighting = StudioLighting()
    property_inspector: PropertyInspector = PropertyInspector()  # NEW
```

**Result:**
- JSONForms in Settings panel sees `enabled_properties` with `x-custom-type="property-inspector"`
- Triggers custom renderer instead of standard array renderer
- Other boolean fields render normally

---

### Phase 2: Frontend Core Infrastructure

#### 2.1 Types & Interfaces

```typescript
// app/src/types/property-inspector.ts

export interface PropertyMetadata {
  dtype: string;
  shape: number[];
  type: "array" | "scalar";
}

export interface PropertyInfo {
  key: string;
  metadata: PropertyMetadata;
  category: "per-particle" | "global";
  prefix: string; // "arrays", "calc", "info", etc.
  enabled: boolean;
}

export interface PropertyValue {
  key: string;
  value: any; // Actual data from API
  formattedValue: string; // For display
  isLoading: boolean;
  error?: string;
}

export interface PropertyInspectorState {
  selectedProperties: Set<string>;
  expandedCategories: Set<string>;
  searchQuery: string;
  sortBy: "name" | "category" | "type";
  filterBy: "all" | "per-particle" | "global";
}
```

#### 2.2 API Client Utilities

```typescript
// app/src/myapi/client.ts

/**
 * Categorize properties based on metadata shape
 */
export const categorizeProperties = (
  metadata: FrameMetadata,
  particleCount: number
): { perParticle: PropertyInfo[]; global: PropertyInfo[] } => {
  const perParticle: PropertyInfo[] = [];
  const global: PropertyInfo[] = [];
  
  metadata.keys.forEach((key) => {
    const meta = metadata.metadata[key];
    if (!meta) return;
    
    // Per-particle detection: first dimension equals particle count
    const isPerParticle = 
      meta.shape.length > 0 && 
      meta.shape[0] === particleCount;
    
    const [prefix] = key.split(".");
    
    const info: PropertyInfo = {
      key,
      metadata: meta,
      category: isPerParticle ? "per-particle" : "global",
      prefix: prefix || "other",
      enabled: false,
    };
    
    if (isPerParticle) {
      perParticle.push(info);
    } else {
      global.push(info);
    }
  });
  
  return { perParticle, global };
};
```

**Note:** Existing `getFrames` API can be used to fetch property values by keys.

#### 2.3 React Query Hooks

```typescript
// app/src/hooks/usePropertyInspector.ts

import { useQuery, useQueries } from "@tanstack/react-query";
import { 
  getFrameMetadata, 
  getPropertyData, 
  categorizeProperties,
  type PropertyInfo,
  type PropertyValue
} from "../myapi/client";

/**
 * Hook to get available properties with categorization
 */
export const useAvailableProperties = (
  roomId: string,
  particleCount: number
) => {
  return useQuery({
    queryKey: ["properties", "available", roomId],
    queryFn: async () => {
      const metadata = await getFrameMetadata(roomId);
      return categorizeProperties(metadata, particleCount);
    },
    enabled: !!roomId && particleCount > 0,
    staleTime: 5 * 60 * 1000, // 5 minutes - metadata rarely changes
  });
};

/**
 * Hook to fetch selected property values
 * Uses parallel queries for optimal performance
 */
export const usePropertyValues = (
  roomId: string,
  frameId: number,
  propertyKeys: string[],
  hoveredParticleId: number | null
) => {
  return useQueries({
    queries: propertyKeys.map((key) => ({
      queryKey: ["property", "value", roomId, frameId, key],
      queryFn: async () => {
        const data = await getPropertyData({
          roomId,
          frameId,
          propertyKeys: [key],
        });
        return {
          key,
          value: data[key]?.value,
          metadata: data[key]?.metadata,
        };
      },
      enabled: !!roomId && propertyKeys.length > 0,
      staleTime: 1000, // Cache for 1 second during rapid frame changes
    })),
  });
};

```

**Note:** No separate settings hook needed - handled by JSONForms in Settings panel.

---

### Phase 3: Custom JSONForms Renderer

This is the **core component** that renders inside the Settings panel when JSONForms encounters `x-custom-type="property-inspector"`.

#### 3.1 PropertyInspectorRenderer (Main Custom Renderer)

**File:** `app/src/components/jsonforms-renderers/PropertyInspectorRenderer.tsx`

This component follows the JSONForms **ControlProps** pattern like `CustomColorPicker` and `DynamicEnumRenderer`.

```typescript
import { withJsonFormsControlProps } from "@jsonforms/react";
import { rankWith, schemaMatches, and, uiTypeIs, ControlProps } from "@jsonforms/core";
import { 
  Box, 
  Typography, 
  Divider, 
  CircularProgress,
  Paper,
  Alert,
} from "@mui/material";
import { useAppStore } from "../../store";
import { useAvailableProperties } from "../../hooks/usePropertyInspector";
import PropertySelector from "./components/PropertySelector";
import PropertyTable from "./components/PropertyTable";

/**
 * PropertyInspectorRenderer - Custom JSONForms control for property inspection
 * 
 * Renders when schema has x-custom-type="property-inspector"
 * Displays property selector + table for per-particle and global properties
 * 
 * ControlProps provided by withJsonFormsControlProps:
 * - data: string[] - Array of selected property keys
 * - handleChange: (path, value) => void - Update form data
 * - path: string - JSON path to this field
 * - label: string - Field label from schema
 * - schema: any - JSON schema for this field
 * - required: boolean - Whether field is required
 * - errors: string - Validation errors
 */
const PropertyInspectorRenderer = ({
  data,
  handleChange,
  path,
  label,
  schema,
  errors,
}: ControlProps) => {
  const { roomId, particleCount, currentFrame, hoveredParticleId } = useAppStore();

  // Fetch available properties with categorization
  const { 
    data: categories, 
    isLoading, 
    isError 
  } = useAvailableProperties(roomId!, particleCount);

  // data is the array of selected property keys (from enabled_properties field)
  const selectedKeys: string[] = data || [];

  // Handler to update JSONForms data
  const handlePropertyToggle = (key: string) => {
    const keys = new Set(selectedKeys);
    if (keys.has(key)) {
      keys.delete(key);
    } else {
      keys.add(key);
    }
    handleChange(path, Array.from(keys));
  };

  const handleClearAll = () => {
    handleChange(path, []);
  };

  if (isLoading) {
    return (
      <Box sx={{ p: 2, display: "flex", justifyContent: "center" }}>
        <CircularProgress size={24} />
        <Typography sx={{ ml: 2 }}>Loading properties...</Typography>
      </Box>
    );
  }

  if (isError) {
    return (
      <Alert severity="error" sx={{ mb: 2 }}>
        Failed to load property metadata
      </Alert>
    );
  }

  return (
    <Paper variant="outlined" sx={{ mb: 2, overflow: "hidden" }}>
      {/* Header */}
      <Box sx={{ p: 2, bgcolor: "action.hover" }}>
        <Typography variant="subtitle2" fontWeight="bold">
          {label || "Property Inspector"}
        </Typography>
        <Typography variant="caption" color="text.secondary">
          Frame {currentFrame} • {particleCount} particles
          {hoveredParticleId !== null && ` • Hovering: Particle ${hoveredParticleId}`}
        </Typography>
      </Box>

      <Divider />

      {/* Property Selector - Collapsible multi-select */}
      <Box sx={{ maxHeight: 300, overflow: "auto" }}>
        <PropertySelector
          perParticleProps={categories?.perParticle || []}
          globalProps={categories?.global || []}
          selectedKeys={selectedKeys}
          onToggle={handlePropertyToggle}
          onClearAll={handleClearAll}
        />
      </Box>

      <Divider />

      {/* Property Table - Shows values for selected properties */}
      <Box sx={{ maxHeight: 400, overflow: "auto" }}>
        <PropertyTable
          propertyKeys={selectedKeys}
          frameId={currentFrame}
          hoveredParticleId={hoveredParticleId}
        />
      </Box>

      {errors && (
        <Alert severity="error" sx={{ m: 1 }}>
          {errors}
        </Alert>
      )}
    </Paper>
  );
};

/**
 * Tester function - determines when to use this renderer
 * Priority 10 to override default array renderers
 */
export const propertyInspectorTester = rankWith(
  10,
  and(
    schemaMatches((schema) => (schema as any)["x-custom-type"] === "property-inspector"),
    uiTypeIs("Control")
  )
);

// Export wrapped with JSONForms HOC
export default withJsonFormsControlProps(PropertyInspectorRenderer);
```

**Key Points:**
- ✅ Follows `ControlProps` pattern like `CustomColorPicker`
- ✅ `data` prop contains current `enabled_properties` array
- ✅ `handleChange(path, value)` updates JSONForms state
- ✅ Renders inside Settings panel automatically
- ✅ Auto-saves via Settings panel's debounced submission

#### 3.2 PropertySelector (Sub-component)

**File:** `app/src/components/jsonforms-renderers/components/PropertySelector.tsx`

Simplified version - just a multi-select list with search:

```typescript
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
  IconButton,
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

  const filterProps = (props: PropertyInfo[]) => {
    if (!searchQuery) return props;
    return props.filter((p) =>
      p.key.toLowerCase().includes(searchQuery.toLowerCase())
    );
  };

  const filteredPerParticle = useMemo(
    () => filterProps(perParticleProps),
    [perParticleProps, searchQuery]
  );

  const filteredGlobal = useMemo(
    () => filterProps(globalProps),
    [globalProps, searchQuery]
  );

  return (
    <Box sx={{ p: 1 }}>
      {/* Search Bar */}
      <TextField
        fullWidth
        size="small"
        placeholder="Search properties..."
        value={searchQuery}
        onChange={(e) => setSearchQuery(e.target.value)}
        InputProps={{
          startAdornment: <SearchIcon sx={{ mr: 1, color: "text.secondary" }} />,
        }}
        sx={{ mb: 1 }}
      />

      {/* Stats Chips */}
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
        {/* Per-Particle Category */}
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
                  secondary={`${prop.metadata.dtype} [${prop.metadata.shape.join(" × ")}]`}
                  primaryTypographyProps={{ fontFamily: "monospace", fontSize: "0.85rem" }}
                />
              </ListItemButton>
            </ListItem>
          ))}
        </Collapse>

        {/* Global Properties Category */}
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
                  secondary={`${prop.metadata.dtype} [${prop.metadata.shape.join(" × ")}]`}
                  primaryTypographyProps={{ fontFamily: "monospace", fontSize: "0.85rem" }}
                />
              </ListItemButton>
            </ListItem>
          ))}
        </Collapse>
      </List>
    </Box>
  );
}
```

#### 3.3 PropertyTable (Simplified MUI-X DataGrid)

**File:** `app/src/components/jsonforms-renderers/components/PropertyTable.tsx`

```typescript
import { useMemo } from "react";
import { Box, Typography } from "@mui/material";
import { DataGrid, GridColDef, GridRowsProp } from "@mui/x-data-grid";
import { usePropertyValues } from "../../../hooks/usePropertyInspector";
import { useAppStore } from "../../../store";

interface PropertyTableProps {
  propertyKeys: string[];
  frameId: number;
  hoveredParticleId: number | null;
}

export default function PropertyTable({
  propertyKeys,
  frameId,
  hoveredParticleId,
}: PropertyTableProps) {
  const { roomId } = useAppStore();
  
  // Fetch property values using React Query
  const propertyQueries = usePropertyValues(
    roomId!,
    frameId,
    propertyKeys,
    hoveredParticleId
  );

  const columns: GridColDef[] = useMemo(
    () => [
      {
        field: "key",
        headerName: "Property",
        flex: 1,
        minWidth: 150,
      },
      {
        field: "value",
        headerName: hoveredParticleId !== null 
          ? `Value (Particle ${hoveredParticleId})` 
          : "Value",
        flex: 2,
        minWidth: 200,
      },
      {
        field: "dtype",
        headerName: "Type",
        width: 80,
      },
      {
        field: "shape",
        headerName: "Shape",
        width: 100,
      },
    ],
    [hoveredParticleId]
  );

  const rows: GridRowsProp = useMemo(() => {
    return propertyQueries.map((query, index) => {
      const key = propertyKeys[index];
      const value = query.data?.value;
      const metadata = query.data?.metadata;

      // Format value for display
      let displayValue = "—";
      if (value !== undefined && value !== null) {
        if (hoveredParticleId !== null && Array.isArray(value)) {
          // Extract particle-specific value
          const particleValue = value[hoveredParticleId];
          displayValue = Array.isArray(particleValue)
            ? `[${particleValue.map(v => typeof v === "number" ? v.toFixed(3) : v).join(", ")}]`
            : String(particleValue);
        } else if (Array.isArray(value)) {
          // Show array preview
          displayValue = value.length <= 5
            ? JSON.stringify(value)
            : `[${value.slice(0, 3).join(", ")}, ... +${value.length - 3} more]`;
        } else {
          displayValue = String(value);
        }
      }

      return {
        id: key,
        key,
        value: query.isLoading ? "Loading..." : displayValue,
        dtype: metadata?.dtype || "...",
        shape: metadata?.shape?.join(" × ") || "...",
      };
    });
  }, [propertyQueries, propertyKeys, hoveredParticleId]);

  if (propertyKeys.length === 0) {
    return (
      <Box sx={{ p: 3, textAlign: "center" }}>
        <Typography color="text.secondary">
          Select properties above to display their values
        </Typography>
      </Box>
    );
  }

  return (
    <DataGrid
      rows={rows}
      columns={columns}
      density="compact"
      disableRowSelectionOnClick
      hideFooter
      autoHeight
      sx={{
        border: "none",
        "& .MuiDataGrid-cell": {
          fontFamily: "monospace",
          fontSize: "0.85rem",
        },
      }}
    />
  );
}
```

**Note:** Simplified - no separate PropertyValueCell component needed for MVP.

#### 3.4 Register Custom Renderer

**Update:** `app/src/utils/jsonforms.ts`

```typescript
import PropertyInspectorRenderer, {
  propertyInspectorTester,
} from "../components/jsonforms-renderers/PropertyInspectorRenderer";

export const customRenderers = [
  ...materialRenderers,
  { tester: dynamicEnumTester, renderer: DynamicEnumRenderer }, // Priority 10
  { tester: propertyInspectorTester, renderer: PropertyInspectorRenderer }, // Priority 10 - NEW
  { tester: customDynamicEnumWithColorPickerTester, renderer: CustomDynamicEnumWithColorPicker }, // Priority 10
  { tester: customColorPickerTester, renderer: CustomColorPicker }, // Priority 5
  { tester: customRangeSliderTester, renderer: CustomRangeSlider },
];
```

**Result:**
- PropertyInspectorRenderer automatically renders in Settings panel
- When JSONForms encounters `enabled_properties` field with `x-custom-type="property-inspector"`
- Fully integrated with Settings auto-save and debouncing

---

### Phase 4: Advanced Features (Future Enhancements)

#### 3.1 Property Comparison Mode
- Select multiple frames
- Display side-by-side comparison
- Highlight differences

#### 3.2 Property Plotting
- Quick plot button for numeric properties
- Histogram for per-particle properties
- Time series for frame evolution

#### 3.3 Property Export
- CSV export with headers
- JSON export for scripting
- Copy to clipboard

#### 3.4 Property Search & Filtering
- Regular expression support
- Type-based filtering (float, int, bool)
- Range filters for numeric values

#### 3.5 Custom Property Expressions
- Compute derived properties: `sqrt(vx² + vy² + vz²)`
- Save as virtual properties
- Apply to selections

---

## User Workflow

### Access Property Inspector

1. **Open Settings panel** (click Settings icon in left sidebar)
2. **Scroll to "Property Inspector" section** (rendered by JSONForms)
3. **Toggle properties on/off** using checkboxes in the property selector
4. **View values in table** below the selector
5. **Hover over particles** in 3D view to see particle-specific values

### Integration with Hover System

- When `hoveredParticleId !== null`:
  - Table header updates: "Value (Particle 42)"
  - Per-particle properties show individual particle values
  - Global properties remain unchanged

---

## Performance Considerations

### 1. **Lazy Loading**
- Only fetch selected properties
- Use React Query's parallel queries
- Cache results per frame

### 2. **Pagination**
- Limit initial property list to 50 items
- Virtual scrolling for large datasets
- Load more on demand

### 3. **Debouncing**
- Debounce search input (300ms)
- Debounce frame changes (100ms)
- Prevent rapid refetches

### 4. **Memoization**
- Memoize property categorization
- Cache formatted values
- Reuse grid columns

### 5. **Smart Caching**
```typescript
// Cache strategy
const queryOptions = {
  staleTime: {
    metadata: 5 * 60 * 1000,    // 5 min - rarely changes
    propertyValue: 1000,         // 1 sec - changes with frame
    propertyList: 10 * 60 * 1000, // 10 min - static per room
  },
  cacheTime: {
    metadata: 30 * 60 * 1000,    // 30 min
    propertyValue: 5 * 60 * 1000, // 5 min
    propertyList: 60 * 60 * 1000, // 1 hour
  },
};
```

---

## Testing Strategy

### Unit Tests
- Property categorization logic
- Value formatting functions
- Filter & search algorithms

### Integration Tests
- Property selection flow
- Data fetching pipeline
- Store updates

### E2E Tests
- Full user workflow
- Keyboard shortcuts
- Error handling

---

## Implementation Phases

### Phase 1: Backend Pydantic Model (Day 1)
- [ ] Add `PropertyInspector` class to `src/zndraw/settings.py`
- [ ] Add to `RoomConfig` class
- [ ] Implement `model_json_schema()` to inject `x-custom-type="property-inspector"`
- [ ] Run `python -m zndraw.settings` to regenerate TypeScript types
- [ ] Test schema generation

### Phase 2: Frontend Types & Hooks (Day 2)
- [ ] Create `app/src/types/property-inspector.ts` with interfaces
- [ ] Add `categorizeProperties()` to `app/src/myapi/client.ts`
- [ ] Create `app/src/hooks/usePropertyInspector.ts` with React Query hooks
  - `useAvailableProperties()` - Fetch and categorize
  - `usePropertyValues()` - Fetch selected property data
- [ ] Test hooks in isolation

### Phase 3: Custom JSONForms Renderer (Days 3-4)
- [ ] Create `app/src/components/jsonforms-renderers/PropertyInspectorRenderer.tsx`
  - Implement tester function
  - Implement renderer component (ControlProps pattern)
  - Handle `data` (string[]) and `handleChange(path, value)`
- [ ] Create `app/src/components/jsonforms-renderers/components/PropertySelector.tsx`
  - Multi-select list with search
  - Collapsible categories (per-particle / global)
- [ ] Create `app/src/components/jsonforms-renderers/components/PropertyTable.tsx`
  - MUI-X DataGrid with 4 columns
  - Handle hover state for particle-specific values
- [ ] Register renderer in `app/src/utils/jsonforms.ts`

### Phase 4: Testing & Polish (Day 5)
- [ ] Test in Settings panel - verify auto-save
- [ ] Test property selection persistence
- [ ] Test hover integration (particle-specific values)
- [ ] Test with large property lists (>100 properties)
- [ ] Add error handling for failed metadata fetch
- [ ] Polish UI spacing and styling

### Phase 5: Documentation & Review (Day 6)
- [ ] Add docstrings to all new functions
- [ ] Update user documentation
- [ ] Create implementation review notes
- [ ] Test edge cases (no properties, all selected, etc.)

---

## Future Considerations

1. **Property History**
   - Track property changes over time
   - Undo/redo for property selection
   - Session persistence

2. **Collaborative Features**
   - Share property selections
   - Broadcast to other users
   - Synchronized inspection

3. **Advanced Visualization**
   - Inline sparklines for numeric properties
   - Color-coded cells by value range
   - Property correlation heatmap

4. **Property Annotations**
   - Add notes to properties
   - Flag important properties
   - Custom property groups

5. **Performance Profiling**
   - Track fetch times
   - Identify slow properties
   - Optimize query patterns

---

## Dependencies

### New Packages Required
```json
{
  "@mui/x-data-grid": "^7.0.0"  // For PropertyTable component
}
```

### Existing Dependencies (Already Available)
- `@tanstack/react-query` - For usePropertyValues hooks
- `@jsonforms/react` - For custom renderer integration
- `@mui/material` - For UI components
- `lodash` - For debouncing (if needed)

### No Breaking Changes
- ✅ Uses existing `getFrames` API
- ✅ Uses existing `getFrameMetadata` API  
- ✅ Integrates with Settings panel (no new sidebar panel)
- ✅ Follows established JSONForms renderer pattern
- ✅ Compatible with current Pydantic/settings architecture

---

## File Structure Summary

```
Backend:
  src/zndraw/settings.py
    + PropertyInspector class
    + Update RoomConfig

Frontend:
  app/src/types/
    + property-inspector.ts (PropertyInfo, PropertyMetadata)
  
  app/src/myapi/client.ts
    + categorizeProperties() function
  
  app/src/hooks/
    + usePropertyInspector.ts (React Query hooks)
  
  app/src/components/jsonforms-renderers/
    + PropertyInspectorRenderer.tsx (main renderer)
    + components/
      + PropertySelector.tsx (sub-component)
      + PropertyTable.tsx (sub-component)
  
  app/src/utils/jsonforms.ts
    + Register propertyInspectorTester/renderer
```

---

## Key Architectural Decisions

### ✅ Why JSONForms Custom Renderer?

**Instead of separate sidebar panel:**
- 🎯 **Consistency:** Follows existing pattern (DynamicEnumRenderer, CustomColorPicker)
- 🎯 **Auto-save:** Leverages Settings panel's debounced submission
- 🎯 **Persistence:** Settings stored per-room via existing API
- 🎯 **Simplicity:** No new navigation items, panel management, or routing

**Pattern Compliance:**
- Uses `withJsonFormsControlProps` HOC
- Uses `rankWith` + `schemaMatches` tester
- Handles `data` (array of strings) and `handleChange(path, value)`
- Renders inside existing Settings panel UI flow

### ✅ Why x-custom-type="property-inspector"?

**Instead of format="property-inspector":**
- `format` is for primitive types (color, range, etc.)
- `x-custom-type` is for complex/composite components
- Matches existing `x-custom-type="dynamic-enum"` pattern

### ✅ Why Categorize Per-Particle vs Global?

**Detection via shape metadata:**
```typescript
isPerParticle = shape[0] === particleCount
```

**Rationale:**
- ✅ Accurate: First dimension = particle count → per-particle property
- ✅ Automatic: No manual tagging required
- ✅ Flexible: Works with any property structure

---

## Conclusion

This implementation provides:
- ✅ **Discoverability** - All properties visible and searchable in Settings
- ✅ **Performance** - Lazy loading via React Query + smart caching
- ✅ **Flexibility** - User-controlled property selection with persistence
- ✅ **Integration** - Seamless fit with JSONForms renderer pattern
- ✅ **Scalability** - Handles hundreds of properties efficiently
- ✅ **Maintainability** - Follows KISS, DRY, SOLID principles
- ✅ **Extensibility** - Easy to add features (export, plotting, filtering)

The Property Inspector becomes a native part of the Settings panel, providing powerful property exploration while maintaining ZnDraw's clean architecture.
