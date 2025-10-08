# Custom JSON Schema Pattern for ZnDraw Geometries

## 🎯 Goal

Design a new JSON schema pattern (and supporting frontend logic) that replaces legacy markers (x-dynamic-enum, x-color-picker, etc.) with:
- A single `x-custom-type` for the conceptual field type
- A flexible `x-features` array for composable UI behaviors
- Renderer logic that handles dynamic visibility (hide if static arrays, etc.)

**Backwards compatibility is not required.**

### 🔄 Migration Status

| Component | Status | Files Affected |
|-----------|--------|----------------|
| Sphere, Bond, Curve | ✅ Migrated | Uses x-custom-type + x-features |
| Arrow | ❌ Needs Migration | Still uses anyOf unions |
| Frontend Renderer | 🚧 Partial | CustomDynamicEnumWithColorPicker exists |
| Unified Renderer | ❌ Not Implemented | Need to create DynamicEnumRenderer |

---

## 📊 Current Implementation Status (API Analysis)

### ✅ Already Migrated to New Pattern

The following geometries **already use** `x-custom-type` and `x-features`:

**Sphere, Bond, Curve:**
```json
{
  "color": {
    "type": "string",
    "x-custom-type": "dynamic-enum",
    "x-features": ["color-picker", "dynamic-atom-props", "free-solo"]
  },
  "position": {
    "type": "string", 
    "x-custom-type": "dynamic-enum",
    "x-features": ["dynamic-atom-props"]
  },
  "radius": {
    "type": "string",
    "x-custom-type": "dynamic-enum",
    "x-features": ["dynamic-atom-props"]  //except for bond.radius which is float
  }
}
```

### ❌ Not Yet Migrated

**Arrow geometry** still uses `anyOf` unions (no x-custom-type markers):
```json
{
  "color": {
    "anyOf": [
      {"type": "string"},
      {"type": "number"},
      {"type": "array", "minItems": 3, "maxItems": 3},
      {"type": "array", "items": {"type": "array"}}
    ],
    "default": "arrays.colors"
  }
}
```

**Key Insight:** Arrow needs to be migrated to match Sphere/Bond/Curve pattern.

🧱 1. Conceptual Model

We’ll base everything on semantic type + feature composition:

Concept	Purpose	Example	Current Usage
x-custom-type	Defines what kind of field it is conceptually	"dynamic-enum"	✅ Sphere, Bond, Curve
x-features	Array of UI behaviors to enable	["dynamic-atom-props", "free-solo", "color-picker"]	✅ Sphere, Bond, Curve

### Feature Definitions

| Feature | Purpose | Behavior |
|---------|---------|----------|
| `dynamic-atom-props` | Populate dropdown from frame metadata | Renderer fetches `metadata.keys` and shows as autocomplete options |
| `free-solo` | Allow custom text input | MUI Autocomplete `freeSolo={true}`, user can type arbitrary values |
| `color-picker` | Add color picker UI | Renders color input alongside autocomplete for hex color selection |

---

## 🧩 2. Schema Design

### Python Implementation Pattern

**Example: Adding to Arrow geometry**

```python
# src/zndraw/geometries/arrow.py

class Arrow(BaseGeometry):
    """Arrow geometry with direction vector."""

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)
        
        # Position field: dropdown of atom props only
        schema["properties"]["position"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["position"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["position"]["type"] = "string"
        schema["properties"]["position"].pop("anyOf", None)
        
        # Color field: dropdown + free text + color picker
        schema["properties"]["color"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["color"]["x-features"] = [
            "color-picker",
            "dynamic-atom-props", 
            "free-solo"
        ]
        schema["properties"]["color"]["type"] = "string"
        schema["properties"]["color"].pop("anyOf", None)
        
        # Radius field: dropdown of atom props (scalars)
        schema["properties"]["radius"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["radius"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["radius"]["type"] = "string"
        schema["properties"]["radius"].pop("anyOf", None)
        
        return schema
```

### JSON Schema Output

```json
{
  "type": "object",
  "properties": {
    "position": {
      "type": "string",
      "description": "Position [x,y,z]. String for dynamic data key, tuple/list for static values.",
      "x-custom-type": "dynamic-enum",
      "x-features": ["dynamic-atom-props"],
      "default": "arrays.positions"
    },
    "color": {
      "type": "string",
      "description": "Color [r,g,b]. String for dynamic data key, tuple/list for static values.",
      "x-custom-type": "dynamic-enum",
      "x-features": ["color-picker", "dynamic-atom-props", "free-solo"],
      "default": "arrays.colors"
    },
    "radius": {
      "type": "string",
      "x-custom-type": "dynamic-enum",
      "x-features": ["dynamic-atom-props"],
      "default": "0.05"
    }
  }
}
```



---

## ⚙️ 3. Frontend Renderer Implementation

### 3.1. Tester Function

```typescript
// app/src/components/jsonforms-renderers/DynamicEnumRenderer.tsx

import { rankWith, schemaMatches, uiTypeIs, and } from "@jsonforms/core";

export const dynamicEnumTester = rankWith(
  10, // High priority
  and(
    schemaMatches((schema) => schema["x-custom-type"] === "dynamic-enum"),
    uiTypeIs("Control")
  )
);
```

### 3.2. Renderer Component Structure

```typescript
// app/src/components/jsonforms-renderers/DynamicEnumRenderer.tsx

import { withJsonFormsControlProps, ControlProps } from "@jsonforms/react";
import { Autocomplete, TextField, Box } from "@mui/material";

const DynamicEnumRenderer = ({
  data,
  handleChange,
  path,
  label,
  schema,
  required,
  errors,
}: ControlProps) => {
  // Extract features
  const features = schema["x-features"] || [];
  const hasColorPicker = features.includes("color-picker");
  const hasFreeSolo = features.includes("free-solo");
  const hasDynamicProps = features.includes("dynamic-atom-props");

  // Get options from injected enum or empty array
  const options = schema.enum || [];
  
  // Detect if current value is static (array/number)
  const isStaticValue = Array.isArray(data) || typeof data === "number";
  
  if (isStaticValue) {
    return <StaticValueDisplay value={data} label={label} />;
  }

  return (
    <Box sx={{ display: "flex", gap: 1, alignItems: "flex-start" }}>
      <Autocomplete
        freeSolo={hasFreeSolo}
        options={options}
        value={data || schema.default || ""}
        onChange={(_, newValue) => handleChange(path, newValue || "")}
        onInputChange={(_, newInputValue, reason) => {
          if (reason === "input" && hasFreeSolo) {
            handleChange(path, newInputValue);
          }
        }}
        renderInput={(params) => (
          <TextField
            {...params}
            label={label}
            required={required}
            error={!!errors}
            helperText={errors}
          />
        )}
        fullWidth
      />
      
      {hasColorPicker && (
        <input
          type="color"
          value={typeof data === "string" && data.startsWith("#") ? data : "#000000"}
          onChange={(e) => handleChange(path, e.target.value)}
          style={{
            width: "50px",
            height: "40px",
            marginTop: "8px",
            border: "1px solid rgba(0, 0, 0, 0.23)",
            borderRadius: "4px",
            cursor: "pointer",
          }}
        />
      )}
    </Box>
  );
};

export default withJsonFormsControlProps(DynamicEnumRenderer);
```

### 3.3. Updated Enum Injection

```typescript
// app/src/utils/jsonforms.ts

/**
 * Inject dynamic enum values based on x-custom-type and x-features.
 */
export const injectDynamicEnums = (
  schema: any,
  metadata: FrameMetadata | undefined
): any => {
  const newSchema = JSON.parse(JSON.stringify(schema));

  const traverse = (obj: any) => {
    if (obj && typeof obj === "object") {
      // Check for x-custom-type="dynamic-enum" with "dynamic-atom-props" feature
      if (
        obj["x-custom-type"] === "dynamic-enum" &&
        Array.isArray(obj["x-features"]) &&
        obj["x-features"].includes("dynamic-atom-props") &&
        metadata?.keys
      ) {
        // Filter keys by expected shape if specified
        let availableKeys = metadata.keys;
        
        obj.enum = availableKeys;
      }

      // Traverse nested objects
      Object.keys(obj).forEach((key) => traverse(obj[key]));
    }
  };

  traverse(newSchema);
  return newSchema;
};
```

### 3.4. Register Renderer

```typescript
// app/src/utils/jsonforms.ts

import DynamicEnumRenderer, {
  dynamicEnumTester,
} from "../components/jsonforms-renderers/DynamicEnumRenderer";

export const customRenderers = [
  ...materialRenderers,
  { tester: dynamicEnumTester, renderer: DynamicEnumRenderer }, // Priority 10
  { tester: customColorPickerTester, renderer: CustomColorPicker }, // Priority 5
  { tester: customRangeSliderTester, renderer: CustomRangeSlider },
];
```

⸻

---

## � Files Requiring Changes

### Python Backend (Schema Generation)

#### ✅ Already Updated (No Changes Needed)
- `src/zndraw/geometries/sphere.py` - Uses x-custom-type and x-features
- `src/zndraw/geometries/bonds.py` - Uses x-custom-type and x-features
- `src/zndraw/geometries/curve.py` - Uses x-custom-type and x-features

#### ❌ Needs Migration
- `src/zndraw/geometries/arrow.py` - **Missing `model_json_schema()` override**
  - Currently relies on Pydantic's default schema generation (anyOf unions)
  - Needs to add x-custom-type and x-features like other geometries

#### 🔍 May Need Review
- `src/zndraw/geometries/base.py` - Base class, check if default behavior can be improved
- `src/zndraw/extensions/analysis.py` - Uses old x-dynamic-enum marker (line 54)

### TypeScript Frontend (Renderer Logic)

#### ❌ Needs Refactoring
- `app/src/utils/jsonforms.ts`
  - `injectDynamicEnums()` - Currently injects enum values for x-dynamic-enum="AVAILABLE_ATOMS_KEYS"
  - **Needs:** Update to work with x-custom-type and x-features
  
- `app/src/components/jsonforms-renderers/CustomDynamicEnumWithColorPicker.tsx`
  - Current tester: Looks for `x-dynamic-enum` AND `x-color-picker`
  - **Needs:** Update tester to look for `x-custom-type="dynamic-enum"` and `"color-picker" in x-features`
  
- `app/src/components/jsonforms-renderers/CustomColorPicker.tsx`
  - Currently handles standalone color picker
  - **May need:** Integration into unified dynamic-enum renderer

#### 📝 Need New Implementation
- Create `CustomDynamicEnumRenderer.tsx` - Single composable renderer that:
  - Checks x-custom-type="dynamic-enum"
  - Reads x-features array to enable/disable:
    - `dynamic-atom-props`: Populate autocomplete from metadata
    - `free-solo`: Allow custom text input
    - `color-picker`: Add color picker UI

#### 🔧 Needs Updates
- `app/src/components/geometry/GeometryForm.tsx` - Uses injectDynamicEnums
- `app/src/components/SecondaryPanel.tsx` - Uses injectDynamicEnums

---

## 🔄 Detailed Implementation Plan

### Phase 1: Migrate Arrow Geometry ✅ (Backend Only)
**Goal:** Make Arrow consistent with Sphere/Bond/Curve

**Tasks:**
1. Add `model_json_schema()` classmethod to `Arrow`
2. Set x-custom-type and x-features for: position, direction, color, radius, scale
3. Test: `curl http://localhost:5000/api/rooms/testroom/geometries/schemas` shows new markers

**Files:**
- `src/zndraw/geometries/arrow.py`

**Acceptance Criteria:**
- Arrow schema matches Sphere pattern with x-custom-type and x-features
- All Arrow fields that accept dynamic data have appropriate features

### Phase 3: Create Unified Dynamic Enum Renderer 🎨
**Goal:** Single composable renderer replaces CustomDynamicEnumWithColorPicker

**Tasks:**
1. Create `app/src/components/jsonforms-renderers/DynamicEnumRenderer.tsx`:
   ```tsx
   - Tester: schema["x-custom-type"] === "dynamic-enum"
   - Render: MUI Autocomplete
   - Feature: "dynamic-atom-props" → populate options from metadata
   - Feature: "free-solo" → freeSolo={true}
   - Feature: "color-picker" → show color input alongside
   ```

2. Update `app/src/utils/jsonforms.ts`:
   - Modify `injectDynamicEnums()` to work with x-custom-type
   - Or eliminate it entirely if renderer reads metadata directly

3. Remove deprecated renderers:
   - `CustomDynamicEnumWithColorPicker.tsx` (replaced by unified renderer)

4. Update renderer registration in `customRenderers` array

**Files:**
- New: `app/src/components/jsonforms-renderers/DynamicEnumRenderer.tsx`
- Update: `app/src/utils/jsonforms.ts`
- Delete: `app/src/components/jsonforms-renderers/CustomDynamicEnumWithColorPicker.tsx`

**Acceptance Criteria:**
- Single renderer handles all x-custom-type="dynamic-enum" fields
- Color picker appears only when "color-picker" in x-features
- Free solo works when "free-solo" in x-features
- Dynamic options populate when "dynamic-atom-props" in x-features

---

### Phase 4: Handle Static vs Dynamic Values 🔀
**Goal:** Hide/disable renderer when value is array/tuple (static data)

**Tasks:**
1. Add value type detection to DynamicEnumRenderer:
   ```tsx
   const isStaticValue = Array.isArray(data) || typeof data === "number"
   if (isStaticValue) {
     return <StaticValueDisplay value={data} />
   }
   ```

---

### Phase 5: Clean Up Legacy Code 🧹
**Goal:** Remove old x-dynamic-enum and x-color-picker markers

**Tasks:**
1. Search and remove any remaining x-dynamic-enum usage
2. Remove x-color-picker standalone usage
3. Update `src/zndraw/extensions/analysis.py` (uses old markers)
4. Update any tests that reference old markers

**Files:**
- `src/zndraw/extensions/analysis.py`
- Search all `*.py` and `*.ts` files for old markers

**Acceptance Criteria:**
- No x-dynamic-enum or x-color-picker in codebase
- All tests pass
- All geometry forms work correctly

---

## 📊 Summary & Next Steps

### ✅ What's Already Done

1. **Backend Schema Pattern Established**
   - Sphere, Bond, and Curve geometries already use x-custom-type and x-features
   - Pattern is proven and working in production

2. **Frontend Infrastructure Exists**
   - `injectDynamicEnums()` function for metadata injection
   - Custom renderer system with JSONForms
   - Metadata API providing keys for dropdowns

### 🚧 What Needs Work

| Priority | Task | Effort | Impact |
|----------|------|--------|--------|
| 🔴 High | Migrate Arrow geometry to new pattern | 30 min | Consistency |
| 🟡 Medium | Create unified DynamicEnumRenderer | 2-3 hours | Better UX |
| 🟢 Low | Clean up old markers in extensions | 1 hour | Code quality |

### 🎯 Recommended Implementation Order

1. **Start with Arrow migration** (Phase 1) - Quick win, establishes consistency
2. **Implement unified renderer** (Phase 3) - Biggest UX improvement
3. **Add shape filtering** (Phase 2) - Nice to have, requires backend changes
4. **Handle static values** (Phase 4) - Polish, improves edge cases
5. **Clean up legacy code** (Phase 5) - Final housekeeping

### 🔍 Key Design Decisions to Confirm

1. **Shape metadata API**: Do we need a new endpoint to expose array shapes, or can we infer from existing data?
2. **Static value editing**: Should array/tuple values be editable in the form, or only via external data sources?
3. **Feature naming**: Keep "color-picker" or rename to "colorpicker"? (Current API uses hyphenated)
4. **Free-solo defaults**: Should all dynamic-enum fields allow free-solo by default?

### 📝 Files Summary

**Backend (Python):**
- ✅ 3 files already migrated: sphere.py, bonds.py, curve.py
- ❌ 1 file needs migration: arrow.py
- 🔍 1 file needs review: analysis.py (extensions)

**Frontend (TypeScript):**
- ❌ 1 file needs major refactor: CustomDynamicEnumWithColorPicker.tsx
- 🆕 1 file needs creation: DynamicEnumRenderer.tsx
- 🔧 2 files need updates: utils/jsonforms.ts, GeometryForm.tsx, SecondaryPanel.tsx

**Total Estimated Effort:** 8-12 hours for complete implementation

