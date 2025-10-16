# Geometry Data Normalization Refactoring Plan

## Problem Statement

The current geometry system has inconsistent data handling:

**Current Issues:**
```python
# Works
box["rotation"] = (0.0, 0.5, 0.0)
box["color"] = "#FF0000"
box["position"] = [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)]

# Doesn't work
box["rotation"] = [(0.0, 0.5, 0.0), (0.0, 0.0, 0.5)]  # Fails validation
box["size"] = [(1.0, 1.0, 1.0), (0.5, 0.5, 0.5)]      # Fails validation
box["position"] = (0.0, 0.0, 0.0)                      # Not allowed
```

## Solution: Normalize to List Format

**All non-dynamic fields stored internally as `List[...]`**

- Input: Accept both `(x,y,z)` and `[(x,y,z), ...]`
- Storage: Always `[(x,y,z)]` or `[(x,y,z), (x,y,z)]`
- Colors: Always hex strings `["#FF0000"]` or `["#FF0000", "#00FF00"]`
- Dynamic refs like `"arrays.positions"` pass through unchanged

---

## Phase 1: Python Backend

### Files Modified
- `src/zndraw/utils.py` - **Convert RGB to hex in `update_colors_and_radii()`**
- `src/zndraw/geometries/base.py`
- `src/zndraw/geometries/box.py`
- `src/zndraw/geometries/plane.py`
- `src/zndraw/geometries/arrow.py`

### 1.1 Update `update_colors_and_radii()` in `utils.py`

**Convert RGB arrays to hex strings:**

```python
def update_colors_and_radii(atoms: ase.Atoms) -> None:
    """Update the colors and radii of the atoms in-place.

    Colors are stored as hex strings instead of RGB arrays.
    """
    if "colors" not in atoms.arrays:
        # Get RGB colors from jmol_colors
        rgb_colors = np.array(
            [
                jmol_colors[atom.number]
                if atom.number < len(jmol_colors)
                else [0.0, 0.0, 0.0]
                for atom in atoms
            ],
            dtype=np.float32,
        )

        # Convert RGB (0-1 range) to hex strings
        hex_colors = []
        for rgb in rgb_colors:
            r, g, b = (np.clip(rgb, 0, 1) * 255).astype(int)
            hex_colors.append(f'#{r:02x}{g:02x}{b:02x}')

        # Store as list of hex strings (not numpy array)
        atoms.set_array("colors", hex_colors)

    if "radii" not in atoms.arrays:
        radii = covalent_radii[atoms.numbers].astype(np.float32)
        atoms.set_array("radii", radii)
```

**Note:** ASE's `set_array()` expects numpy arrays, so we need to check if it accepts list of strings. If not, we'll store as object array:
```python
atoms.set_array("colors", np.array(hex_colors, dtype=object))
```

### 1.2 Update Type Aliases in `base.py`

Use 1D, 2D, 3D type aliases for clarity:
```python
ColorProp = Union[str, list[str]]  # Always hex strings
SizeProp = Union[str, list[float]]  # For 1D (radius)
Size2DProp = Union[str, list[tuple[float, float]]]  # For Plane
Size3DProp = Union[str, list[tuple[float, float, float]]]  # For Box
```

### 1.3 Add Field Validators

Add to **Box, Plane, Arrow** classes:

```python
@field_validator("position", "rotation", "size", mode="before")
@classmethod
def normalize_vector_fields(cls, v):
    """Normalize vector fields to list of tuples."""
    if v is None:
        return []
    if isinstance(v, str):  # Dynamic reference
        return v
    if isinstance(v, tuple):  # Single tuple -> wrap in list
        return [v]
    return v  # Already a list

@field_validator("color", mode="before")
@classmethod
def normalize_color(cls, v):
    """Normalize color to list of hex strings."""
    if v is None:
        return []
    # Single hex color (not dynamic) -> wrap in list
    if isinstance(v, str) and v.startswith("#"):
        return [v]
    # Dynamic reference -> pass through
    return v
```

### 1.4 Update Defaults

```python
position: ... = Field(default=[(0.0, 0.0, 0.0)], ...)
rotation: ... = Field(default=[(0.0, 0.0, 0.0)], ...)
size: ... = Field(default=[(1.0, 1.0, 1.0)], ...)  # Box
size: ... = Field(default=[(1.0, 1.0)], ...)       # Plane
color: ... = Field(default=["#808080"], ...)
```

---

## Phase 2: TypeScript Data Processing

### 2.1 Remove Color Conversion

**Delete:**
- `hexToRgb()` function from `app/src/utils/colorUtils.ts`

**Keep:**
- `isHexColor()` - still needed for validation
- `shouldFetchAsFrameData()` - still needed for query logic

### 2.2 Simplify Processing Functions in `geometryData.ts`

Replace specialized functions with generic dimensional handlers:

#### processColorData
```typescript
/**
 * Process color data - returns array of hex strings
 * Backend sends: string (dynamic ref) | string[] (hex list)
 * Fetched data: string[] (hex list from server)
 */
export function processColorData(
  propValue: string | string[],
  fetchedValue: any,
  count: number
): string[] {
  // Fetched data is already hex strings from backend
  if (fetchedValue) {
    return fetchedValue as string[];
  }

  // Backend always sends list of hex strings ["#FF0000", ...]
  if (Array.isArray(propValue)) {
    return propValue;
  }

  // Shouldn't reach here - dynamic refs should be fetched
  throw new Error(`Dynamic color reference not fetched: ${propValue}`);
}
```

#### process3DData
```typescript
/**
 * Process 3D vector data (position, rotation, size)
 * Backend sends: string (dynamic ref) | [[x,y,z], [x,y,z], ...]
 * Returns: flat array [x,y,z,x,y,z,...]
 */
export function process3DData(
  propValue: string | number[][],
  fetchedValue: any
): number[] {
  if (fetchedValue) {
    return Array.from(fetchedValue);
  }

  if (Array.isArray(propValue) && Array.isArray(propValue[0])) {
    return propValue.flat();
  }

  return [];
}
```

#### process2DData
```typescript
/**
 * Process 2D vector data (plane size)
 * Backend sends: string (dynamic ref) | [[w,h], [w,h], ...]
 * Returns: flat array [w,h,w,h,...]
 */
export function process2DData(
  propValue: string | number[][],
  fetchedValue: any
): number[] {
  if (fetchedValue) {
    return Array.from(fetchedValue);
  }

  if (Array.isArray(propValue) && Array.isArray(propValue[0])) {
    return propValue.flat();
  }

  return [];
}
```

#### process1DData
```typescript
/**
 * Process 1D scalar data (radius, scale)
 * Backend sends: string (dynamic ref) | number | [v1, v2, ...]
 * Returns: flat array [v1, v2, v3, ...]
 */
export function process1DData(
  propValue: string | number | number[],
  fetchedValue: any,
  count: number
): number[] {
  if (fetchedValue) {
    return Array.from(fetchedValue);
  }

  if (typeof propValue === 'number') {
    return Array(count).fill(propValue);
  }

  if (Array.isArray(propValue)) {
    return propValue;
  }

  return [];
}
```

#### getInstanceCount
```typescript
export function getInstanceCount(
  positionProp: string | number[][],
  fetchedPosition: any
): number {
  if (fetchedPosition) {
    return fetchedPosition.length / 3;
  }

  if (Array.isArray(positionProp)) {
    return positionProp.length;  // Backend always sends [[x,y,z], ...]
  }

  return 0;  // Dynamic reference not yet fetched
}
```

#### getColorCount
```typescript
/**
 * Get color count to detect single-color mode for UI
 */
export function getColorCount(
  colorProp: string | string[]
): number {
  if (typeof colorProp === 'string') {
    return 0;  // Dynamic reference
  }
  return colorProp.length;  // List of hex strings
}
```

### 2.3 Update Box.tsx Color Handling

```typescript
// Process colors
const fetchedColor = typeof colorProp === 'string' ? colorData?.[colorProp as string] : undefined;
const colorHexArray = processColorData(colorProp, fetchedColor, finalCount);

// Handle shared color (single color for all instances)
const colorsAreShared = colorHexArray.length === 1 && finalCount > 1;
const finalColorHex = colorsAreShared
  ? Array(finalCount).fill(colorHexArray[0])
  : colorHexArray;

// Validation - colors should match instance count
const isDataValid = validateArrayLengths(
  { positions: finalPositions, sizes: finalSizes, rotations: finalRotations },
  { positions: finalCount * 3, sizes: finalCount * 3, rotations: finalCount * 3 }
) && (finalColorHex.length === finalCount);

// Set colors using THREE.Color
for (let i = 0; i < finalCount; i++) {
  // ... position, rotation, size logic ...

  // Set color directly from hex string
  _color.set(finalColorHex[i]);
  mainMesh.setColorAt(i, _color);
}
```

---

## Phase 3: Array Editor

### 3.1 Update `arrayEditor.ts`

#### Add Color Field Type
```typescript
export type ArrayFieldType =
  | 'position'
  | 'direction'
  | 'color'
  | 'radius'
  | 'scale'
  | 'rotation'
  | 'size'
  | 'generic';
```

#### Update Field Type Config
```typescript
color: {
  dimensions: 1,  // One hex string per row
  columnLabels: ['Color'],
  supportsSingleValue: true,
  defaultValue: '#808080',  // Hex string
  isColor: true,  // NEW: flag for special color handling
},
rotation: {
  dimensions: 3,
  columnLabels: ['X', 'Y', 'Z'],
  supportsSingleValue: true,  // CHANGED
  defaultValue: 0,
},
size: {
  dimensions: 3,
  columnLabels: ['W', 'H', 'D'],
  supportsSingleValue: true,
  defaultValue: 1.0,
},
```

#### Update normalizeToArray for Colors
```typescript
export function normalizeToArray(
  value: string | number | string[] | number[] | number[][],
  fieldType: ArrayFieldType
): (number | string)[][] {
  const config = getFieldTypeConfig(fieldType);

  // Special handling for colors
  if (fieldType === 'color') {
    if (typeof value === 'string' && value.startsWith('#')) {
      return [[value]];  // Single hex string
    }
    if (Array.isArray(value) && typeof value[0] === 'string') {
      return value.map(hex => [hex]);  // List of hex strings
    }
  }

  // Already 2D array
  if (Array.isArray(value) && Array.isArray(value[0])) {
    return value;
  }

  // 1D array
  if (Array.isArray(value)) {
    if (value.length === config.dimensions) {
      return [value];  // Single row
    }
    // Chunk by dimensions
    const result: any[][] = [];
    for (let i = 0; i < value.length; i += config.dimensions) {
      result.push(value.slice(i, i + config.dimensions));
    }
    return result;
  }

  // Single number
  if (typeof value === 'number') {
    return [Array(config.dimensions).fill(value)];
  }

  // Default
  return [Array(config.dimensions).fill(config.defaultValue)];
}
```

#### Update denormalizeFromArray
```typescript
export function denormalizeFromArray(
  arrayValue: (number | string)[][],
  fieldType: ArrayFieldType
): number | string | number[] | string[] | number[][] {
  const config = getFieldTypeConfig(fieldType);

  if (arrayValue.length === 0) {
    return [];
  }

  // Special handling for colors - always return list of hex strings
  if (fieldType === 'color') {
    return arrayValue.map(row => row[0] as string);
  }

  // For 1D fields (radius, scale)
  if (config.dimensions === 1) {
    if (arrayValue.length === 1 && config.supportsSingleValue) {
      return arrayValue[0][0] as number;
    }
    return arrayValue.map(row => row[0] as number);
  }

  // For multi-dimensional fields - always return list of tuples
  return arrayValue as number[][];
}
```

### 3.2 Update `ArrayEditorDialog.tsx`

#### Color Picker Integration

```typescript
// Add color picker state
const [colorPickerMode, setColorPickerMode] = useState(false);

// DataGrid columns - special handling for colors
const columns: GridColDef<RowData>[] = useMemo(() => {
  const cols: GridColDef<RowData>[] = [
    { field: 'id', headerName: '#', width: 60, editable: false },
  ];

  if (fieldType === 'color') {
    // Color field: show color picker + hex input
    cols.push({
      field: 'col0',
      headerName: 'Color',
      width: 200,
      editable: true,
      renderCell: (params) => (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <input
            type="color"
            value={params.value}
            onChange={(e) => {
              const newRow = { ...params.row, col0: e.target.value };
              processRowUpdate(newRow);
            }}
            style={{ width: 40, height: 30, cursor: 'pointer' }}
          />
          <span style={{ fontFamily: 'monospace' }}>{params.value}</span>
        </Box>
      ),
      renderEditCell: (params) => (
        <TextField
          value={params.value}
          onChange={(e) => params.api.setEditCellValue({
            id: params.id,
            field: 'col0',
            value: e.target.value
          })}
          placeholder="#FF0000"
          fullWidth
        />
      ),
    });
  } else {
    // Numeric fields
    for (let i = 0; i < config.dimensions; i++) {
      cols.push({
        field: `col${i}`,
        headerName: config.columnLabels[i],
        width: 120,
        editable: true,
        type: 'number',
      });
    }
  }

  return cols;
}, [config, fieldType]);
```

#### Validation for Colors
```typescript
export function validateArrayData(
  arrayValue: (number | string)[][],
  fieldType: ArrayFieldType
): { valid: boolean; errors: string[] } {
  const errors: string[] = [];
  const config = getFieldTypeConfig(fieldType);

  arrayValue.forEach((row, idx) => {
    if (row.length !== config.dimensions) {
      errors.push(`Row ${idx + 1} has ${row.length} values, expected ${config.dimensions}`);
    }

    // Color validation
    if (fieldType === 'color') {
      const hex = row[0] as string;
      if (!isHexColor(hex)) {
        errors.push(`Row ${idx + 1}: invalid hex color "${hex}"`);
      }
    }

    // Numeric range validation
    if (config.valueRange && typeof row[0] === 'number') {
      row.forEach((val, colIdx) => {
        if (typeof val === 'number' && (val < config.valueRange![0] || val > config.valueRange![1])) {
          errors.push(
            `Row ${idx + 1}, ${config.columnLabels[colIdx]}: value ${val} outside range`
          );
        }
      });
    }
  });

  return { valid: errors.length === 0, errors };
}
```

---

## Phase 4: UI Components

### 4.1 Update `DynamicEnumRenderer.tsx`

#### Detect Static Color Arrays
```typescript
// Detect if value is static color array
const isStaticValue = Array.isArray(data) || typeof data === "number";
const isColorArray = Array.isArray(data) && data.length > 0 && typeof data[0] === 'string' && data[0].startsWith('#');
```

#### Show Color Picker for Single-Color Arrays
```typescript
// For single-color arrays, show color picker alongside autocomplete
const showColorPicker = hasColorPicker && (
  (typeof data === 'string' && data.startsWith('#')) ||
  (isColorArray && data.length === 1)
);

const colorValue = useMemo(() => {
  if (typeof data === 'string' && data.startsWith('#')) return data;
  if (isColorArray && data.length === 1) return data[0];
  return "#000000";
}, [data, isColorArray]);

const handleColorPickerChange = (hex: string) => {
  // Maintain array format for arrays
  if (Array.isArray(data)) {
    handleChange(path, [hex]);
  } else {
    handleChange(path, hex);
  }
};
```

### 4.2 Update Three.js Components

**Update type definitions:**
```typescript
interface BoxData {
  position: string | number[][];
  size: string | number[][];
  color: string | string[];  // Dynamic ref or list of hex
  rotation: string | number[][];
  material: string;
  scale: number;
  opacity: number;
}
```

**Remove RGB processing, use THREE.Color directly:**
```typescript
// Process colors to hex array
const colorHexArray = processColorData(colorProp, fetchedColor, finalCount);

// Handle shared color
const finalColorHex = colorHexArray.length === 1 && finalCount > 1
  ? Array(finalCount).fill(colorHexArray[0])
  : colorHexArray;

// Set colors
for (let i = 0; i < finalCount; i++) {
  _color.set(finalColorHex[i]);  // THREE.Color.set() accepts hex directly
  mainMesh.setColorAt(i, _color);
}
```

---

## Phase 5: Testing

### 5.1 Python Backend Tests

**File:** `tests/test_geometries.py`

```python
def test_box_normalization():
    """Test Box normalizes inputs to list format"""
    # Single tuple/hex inputs -> wrapped in lists
    box1 = Box(position=(0, 0, 0), rotation=(0, 0.5, 0), color="#FF0000")
    assert box1.position == [(0.0, 0.0, 0.0)]
    assert box1.rotation == [(0.0, 0.5, 0.0)]
    assert box1.color == ["#FF0000"]

    # List inputs -> kept as lists
    box2 = Box(
        position=[(0, 0, 0), (1, 1, 1)],
        rotation=[(0, 0.5, 0), (0, 0, 0.5)],
        color=["#FF0000", "#00FF00"]
    )
    assert len(box2.position) == 2
    assert len(box2.rotation) == 2
    assert box2.color == ["#FF0000", "#00FF00"]

    # Dynamic references -> pass through
    box3 = Box(position="arrays.positions", color="arrays.colors")
    assert box3.position == "arrays.positions"
    assert box3.color == "arrays.colors"

def test_color_conversion():
    """Test that atoms colors are converted to hex"""
    atoms = ase.Atoms('H2O', positions=[(0,0,0), (1,0,0), (0,1,0)])
    update_colors_and_radii(atoms)

    colors = atoms.arrays["colors"]
    assert len(colors) == 3
    assert all(isinstance(c, str) and c.startswith('#') for c in colors)
```

### 5.2 Frontend Tests

**Manual Testing Checklist:**

- [ ] Box with `position=[(0,0,0)]` renders correctly
- [ ] Box with `position=[(0,0,0), (1,1,1)]` renders two boxes
- [ ] Box with `rotation=[(0,0.5,0)]` applies rotation to all instances
- [ ] Box with `rotation=[(0,0.5,0), (0,0,0.5)]` applies per-instance rotations
- [ ] Box with `color=["#FF0000"]` shows single red box with color picker enabled
- [ ] Box with `color=["#FF0000", "#00FF00"]` shows red and green boxes
- [ ] Loading atoms with `arrays.colors` dynamic reference shows correct colors
- [ ] Array editor for rotation shows correct rows
- [ ] Array editor for color shows color pickers per row
- [ ] Editing colors in array editor with hex input works
- [ ] Editing colors with color picker works
- [ ] Saving from array editor preserves list format
- [ ] Dynamic refs `position="arrays.positions"` still work

---

## Implementation Order

1. **Phase 1:** Python changes (2-3h)
   - Convert RGB to hex in `update_colors_and_radii`
   - Add validators to geometry classes
2. **Phase 2:** TypeScript data processing (2-3h)
   - Remove `hexToRgb`
   - Simplify to generic processing functions
   - Update Three.js color handling
3. **Phase 3:** Array editor (2-3h)
   - Add color support
   - Update denormalization
4. **Phase 4:** UI components (2-3h)
   - DynamicEnumRenderer color picker
   - Three.js component updates
5. **Phase 5:** Testing (2-3h)

**Total:** ~10-14 hours

---

## Files Affected

### Python Backend
- `src/zndraw/utils.py` - Convert RGB to hex
- `src/zndraw/geometries/base.py` - Type aliases
- `src/zndraw/geometries/box.py` - Validators, defaults
- `src/zndraw/geometries/plane.py` - Validators, defaults
- `src/zndraw/geometries/arrow.py` - Validators, defaults

### TypeScript Frontend
- `app/src/utils/colorUtils.ts` - Remove `hexToRgb`
- `app/src/utils/geometryData.ts` - Generic processing functions
- `app/src/utils/arrayEditor.ts` - Color support, denormalization
- `app/src/components/jsonforms-renderers/ArrayEditorDialog.tsx` - Color picker rows
- `app/src/components/jsonforms-renderers/DynamicEnumRenderer.tsx` - Color picker for single-element arrays
- `app/src/components/three/Box.tsx` - THREE.Color usage
- `app/src/components/three/Plane.tsx` - THREE.Color usage
- `app/src/components/three/Arrow.tsx` - THREE.Color usage

### Tests
- `tests/test_geometries.py` - Normalization tests

---

## Success Criteria

- [ ] All geometry fields accept both single values and lists
- [ ] All non-dynamic fields stored as lists internally
- [ ] Colors converted to hex in `update_colors_and_radii()`
- [ ] Colors always hex strings throughout pipeline
- [ ] THREE.Color.set() used directly with hex strings
- [ ] Array editor supports color editing with pickers
- [ ] Color picker shows for single-element color arrays
- [ ] No validation errors for list format fields
- [ ] Dynamic references still work
- [ ] All tests pass
