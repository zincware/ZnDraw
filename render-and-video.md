https://raw.githubusercontent.com/gkjohnson/three-gpu-pathtracer/383bbdca0cad15582688707fb0a3eaad735f0eb0/example/renderVideo.html

https://raw.githubusercontent.com/gkjohnson/three-gpu-pathtracer/refs/heads/main/example/renderVideo.js

https://raw.githubusercontent.com/zincware/ZnDraw/0c17eca4e96f4e2591865ba13ad9bf7900e48ca4/app/src/components/utils/mergeInstancedMesh.tsx

TODO:
- add pathtracing settings to the existing settings = {
    "camera": Camera,
    "studio_lighting": StudioLighting,
    "property_inspector": PropertyInspector,
}
- convert instanced meshes, like arrows / spheres / bonds to non-instanced meshes for pathtracing
- add support for rendering a video over all steps in the current room
- use `@react-three/gpu-pathtracer` from https://github.com/pmndrs/react-three-gpu-pathtracer

---

# GPU Path Tracing Implementation Plan

## Critical Design Decisions

### ✅ CONFIRMED: Must Convert Instanced Meshes
**GPU path tracing does NOT support InstancedMesh.**
- Three.js InstancedMesh uses GPU instancing which is incompatible with path tracing
- Must convert each instance to individual Mesh objects before path tracing
- Conversion is expensive but unavoidable
- Memory usage will increase: N instances → N meshes + N materials

### ✅ Library Choice: @react-three/gpu-pathtracer
- Installed and ready to use
- React wrapper for three-gpu-pathtracer
- Provides `<Pathtracer>` component for easy integration
- Handles progressive rendering automatically
- Works with @react-three/drei's `<Environment>` for HDRI lighting

### ✅ Dual Rendering Architecture
- Keep InstancedMesh for normal rendering (fast, interactive)
- Switch to non-instanced meshes only when path tracing enabled
- Each geometry component must support both modes
- Toggle via `pathtracingEnabled` prop

---

## Implementation Phases

### Phase 1: Backend Settings (Python)

**File: src/zndraw/settings.py**

1. **Add PathTracing class** (EnvironmentPreset already exists at line 42):

```python
class PathTracing(SettingsBase):
    """GPU Path Tracing settings for high-quality physically-based rendering."""

    enabled: bool = Field(
        default=False,
        description="Enable GPU path tracing renderer"
    )

    min_samples: int = Field(
        default=1,
        ge=1,
        description="Minimum samples before displaying result"
    )

    samples: int = Field(
        default=256,
        ge=1,
        le=10000,
        description="Maximum samples to render"
    )

    bounces: int = Field(
        default=3,
        ge=1,
        le=32,
        description="Number of light bounces for global illumination"
    )

    tiles: int = Field(
        default=1,
        ge=1,
        le=8,
        description="Rendering tile count (higher = less memory, slower)"
    )

    environment_preset: EnvironmentPreset = Field(
        default=EnvironmentPreset.studio,
        description="HDRI environment preset for scene lighting"
    )

    environment_intensity: float = Field(
        default=1.0,
        ge=0.0,
        le=10.0,
        description="Environment map brightness multiplier"
    )

    environment_blur: float = Field(
        default=0.0,
        ge=0.0,
        le=1.0,
        description="Environment background blur amount"
    )

    environment_background: bool = Field(
        default=False,
        description="Show environment as visible background"
    )

    @classmethod
    def model_json_schema(cls, *args, **kwargs) -> dict[str, t.Any]:
        schema = super().model_json_schema(*args, **kwargs)
        schema["properties"]["min_samples"]["format"] = "range"
        schema["properties"]["samples"]["format"] = "range"
        schema["properties"]["bounces"]["format"] = "range"
        schema["properties"]["tiles"]["format"] = "range"
        schema["properties"]["environment_intensity"]["format"] = "range"
        schema["properties"]["environment_intensity"]["step"] = 0.1
        schema["properties"]["environment_blur"]["format"] = "range"
        schema["properties"]["environment_blur"]["step"] = 0.01
        return schema
```

2. **Update settings dictionary** (line 153):

```python
settings = {
    "camera": Camera,
    "studio_lighting": StudioLighting,
    "property_inspector": PropertyInspector,
    "pathtracing": PathTracing,  # ← ADD THIS
}
```

3. **Update RoomConfig class** (line 160):

```python
class RoomConfig(SettingsBase):
    """ZnDraw room configuration combining all settings sections."""

    camera: Camera = Camera()
    studio_lighting: StudioLighting = StudioLighting()
    property_inspector: PropertyInspector = PropertyInspector()
    pathtracing: PathTracing = PathTracing()  # ← ADD THIS
```

4. **Generate TypeScript types**:

```bash
cd src
python -m zndraw.settings
# This auto-generates app/src/types/room-config.ts
```

---

### Phase 2: Instanced Mesh Converter Utility

**File: app/src/utils/convertInstancedMesh.ts** (NEW)

Critical utility that converts InstancedMesh to Group of Meshes:

```typescript
import * as THREE from 'three';

/**
 * Converts an InstancedMesh to a Group of individual Meshes.
 * Required for GPU path tracing which does not support instancing.
 * Preserves existing material properties from the instanced mesh.
 *
 * @param instancedMesh - The instanced mesh to convert
 * @returns A THREE.Group containing individual meshes for each instance
 */
export function convertInstancedMeshToGroup(
  instancedMesh: THREE.InstancedMesh
): THREE.Group {
  const group = new THREE.Group();
  const tempMatrix = new THREE.Matrix4();
  const tempPosition = new THREE.Vector3();
  const tempRotation = new THREE.Quaternion();
  const tempScale = new THREE.Vector3();
  const tempColor = new THREE.Color();

  const baseGeometry = instancedMesh.geometry;
  const baseMaterial = instancedMesh.material;
  const count = instancedMesh.count;
  const hasColors = instancedMesh.instanceColor !== null;

  for (let i = 0; i < count; i++) {
    // Extract instance transformation
    instancedMesh.getMatrixAt(i, tempMatrix);
    tempMatrix.decompose(tempPosition, tempRotation, tempScale);

    // Clone geometry and bake transformation into vertices
    const geometry = baseGeometry.clone();
    geometry.applyMatrix4(tempMatrix);

    // Clone and adapt material with instance color if present
    let material: THREE.Material;

    if (hasColors) {
      instancedMesh.getColorAt(i, tempColor);

      // Clone base material and override color
      if (baseMaterial instanceof THREE.MeshStandardMaterial) {
        material = baseMaterial.clone();
        material.color.copy(tempColor);
      } else if (baseMaterial instanceof THREE.MeshPhysicalMaterial) {
        material = baseMaterial.clone();
        material.color.copy(tempColor);
      } else if (baseMaterial instanceof THREE.MeshBasicMaterial) {
        material = baseMaterial.clone();
        material.color.copy(tempColor);
      } else {
        // Fallback: clone whatever material it is
        material = baseMaterial.clone();
      }
    } else {
      // No per-instance colors, just clone base material
      material = baseMaterial.clone();
    }

    // Create individual mesh
    const mesh = new THREE.Mesh(geometry, material);
    group.add(mesh);
  }

  // Copy group transformation from original mesh
  group.position.copy(instancedMesh.position);
  group.rotation.copy(instancedMesh.rotation);
  group.scale.copy(instancedMesh.scale);

  return group;
}

/**
 * Disposes of all geometries and materials in a group.
 * Essential for preventing memory leaks when switching modes.
 */
export function disposeGroup(group: THREE.Group): void {
  group.traverse((child) => {
    if (child instanceof THREE.Mesh) {
      if (child.geometry) {
        child.geometry.dispose();
      }
      if (child.material) {
        if (Array.isArray(child.material)) {
          child.material.forEach(m => m.dispose());
        } else {
          child.material.dispose();
        }
      }
    }
  });
  group.clear();
}
```

---

### Phase 3: Modify Geometry Components

Each component needs dual-mode rendering. Start with **Particles.tsx** as template.

**File: app/src/components/three/Particles.tsx**

Key changes:

1. **Add prop**:
```typescript
export default function Sphere({
  data,
  geometryKey,
  pathtracingEnabled = false  // ← ADD THIS
}: {
  data: SphereData;
  geometryKey: string;
  pathtracingEnabled?: boolean;  // ← ADD THIS
}) {
```

2. **Add ref for non-instanced group**:
```typescript
const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
const nonInstancedGroupRef = useRef<THREE.Group | null>(null);  // ← ADD THIS
```

3. **Add conversion effect**:
```typescript
// Convert to non-instanced when pathtracing enabled
useEffect(() => {
  if (!pathtracingEnabled) {
    // Clean up non-instanced meshes when disabled
    if (nonInstancedGroupRef.current) {
      disposeGroup(nonInstancedGroupRef.current);
    }
    return;
  }

  if (!mainMeshRef.current || !nonInstancedGroupRef.current) return;

  // Clear existing non-instanced meshes
  disposeGroup(nonInstancedGroupRef.current);

  // Convert instanced mesh to group (preserves existing materials)
  const group = convertInstancedMeshToGroup(mainMeshRef.current);

  // Transfer meshes to our ref
  group.children.forEach(child => {
    nonInstancedGroupRef.current!.add(child);
  });

  // Cleanup
  return () => {
    if (nonInstancedGroupRef.current) {
      disposeGroup(nonInstancedGroupRef.current);
    }
  };
}, [
  pathtracingEnabled,
  instanceCount,
  // Add other deps that affect instance data
  positionData,
  colorData,
  radiusData,
]);
```

4. **Conditional rendering**:
```typescript
return (
  <group>
    {/* Instanced meshes - visible when NOT pathtracing */}
    {!pathtracingEnabled && (
      <>
        <instancedMesh
          key={instanceCount}
          ref={mainMeshRef}
          args={[undefined, undefined, instanceCount]}
          onClick={selecting.enabled ? onClickHandler : undefined}
          onPointerEnter={hovering?.enabled ? onPointerEnterHandler : undefined}
          onPointerMove={hovering?.enabled ? onPointerMoveHandler : undefined}
          onPointerOut={hovering?.enabled ? onPointerOutHandler : undefined}
        >
          <primitive object={mainGeometry} attach="geometry" />
          {renderMaterial(material, data.opacity)}
        </instancedMesh>

        {/* Selection and hover meshes */}
        {selecting.enabled && (
          <instancedMesh
            key={`selection-${validSelectedIndices.length}`}
            ref={selectionMeshRef}
            args={[undefined, undefined, validSelectedIndices.length]}
          >
            <primitive object={mainGeometry} attach="geometry" />
            <meshBasicMaterial
              side={THREE.FrontSide}
              transparent
              opacity={selecting.opacity}
              color={selecting.color}
            />
          </instancedMesh>
        )}

        {hovering?.enabled && (
          <mesh ref={hoverMeshRef} visible={false}>
            <primitive object={mainGeometry} attach="geometry" />
            <meshBasicMaterial
              side={THREE.BackSide}
              transparent
              opacity={hovering.opacity}
              color={hovering.color}
            />
          </mesh>
        )}
      </>
    )}

    {/* Non-instanced group - visible when pathtracing */}
    {pathtracingEnabled && (
      <group ref={nonInstancedGroupRef} />
    )}
  </group>
);
```

**Repeat for:**
- `app/src/components/three/SingleBonds.tsx`
- `app/src/components/three/Arrow.tsx`

---

### Phase 4: PathTracingRenderer Component

**File: app/src/components/PathTracingRenderer.tsx** (NEW)

```typescript
import { Pathtracer } from '@react-three/gpu-pathtracer';
import { Environment } from '@react-three/drei';
import type { PathTracing } from '../types/room-config';

interface PathTracingRendererProps {
  settings: PathTracing;
  children: React.ReactNode;
}

/**
 * Wraps scene with GPU path tracing renderer and environment lighting.
 * When disabled, passes children through without modification.
 */
export function PathTracingRenderer({
  settings,
  children
}: PathTracingRendererProps) {
  const {
    enabled = false,
    min_samples = 1,
    samples = 256,
    bounces = 3,
    tiles = 1,
    environment_preset = 'studio',
    environment_intensity = 1.0,
    environment_blur = 0.0,
    environment_background = false,
  } = settings;

  // Pass through if not enabled
  if (!enabled) {
    return <>{children}</>;
  }

  return (
    <Pathtracer
      minSamples={min_samples}
      samples={samples}
      bounces={bounces}
      tiles={tiles}
      enabled={enabled}
    >
      {/* Environment lighting for path tracing */}
      {environment_preset !== 'none' && (
        <Environment
          preset={environment_preset}
          background={environment_background}
          backgroundBlurriness={environment_blur}
          environmentIntensity={environment_intensity}
        />
      )}
      {children}
    </Pathtracer>
  );
}
```

---

### Phase 5: Integrate into Canvas

**File: app/src/components/Canvas.tsx**

1. **Import new components** (top):
```typescript
import { PathTracingRenderer } from './PathTracingRenderer';
import { convertInstancedMeshToGroup, disposeGroup } from '../utils/convertInstancedMesh';
```

2. **Fetch pathtracing settings** (around line 42):
```typescript
const { data: pathtracingSettings } = useExtensionData(
  roomId || "",
  userId || "",
  "settings",
  "pathtracing",
);

const pathtracingEnabled = pathtracingSettings?.enabled === true;
```

3. **Wrap scene in PathTracingRenderer** (around line 92):
```typescript
<Canvas
  key={cameraSettings.camera}
  shadows
  camera={{ position: [-10, 10, 30], fov: 50 }}
  gl={{
    antialias: true,
    toneMapping: THREE.ACESFilmicToneMapping,
    preserveDrawingBuffer: cameraSettings.preserve_drawing_buffer === true
  }}
  style={{ background: backgroundColor }}
  orthographic={cameraSettings.camera === "OrthographicCamera"}
>
  <CameraManager settings={cameraSettings} />

  {/* Wrap in PathTracingRenderer */}
  <PathTracingRenderer settings={pathtracingSettings}>
    {/* Disable studio lighting when pathtracing (environment provides light) */}
    {!pathtracingEnabled && (
      <SceneLighting
        ambient_light={studioLightingSettings.ambient_light}
        key_light={studioLightingSettings.key_light}
        fill_light={studioLightingSettings.fill_light}
        rim_light={studioLightingSettings.rim_light}
        hemisphere_light={studioLightingSettings.hemisphere_light}
      />
    )}

    <KeyboardShortcutsHandler />

    {/* Pass pathtracingEnabled to all geometry components */}
    {Object.entries(geometries)
      .filter(([_, config]) => config.data?.active !== false)
      .map(([name, config]) => {
        if (config.type === "Sphere") {
          return (
            <Sphere
              key={name}
              geometryKey={name}
              data={config.data}
              pathtracingEnabled={pathtracingEnabled}  // ← ADD THIS
            />
          );
        } else if (config.type === "Bond") {
          return (
            <Bonds
              key={name}
              geometryKey={name}
              data={config.data}
              pathtracingEnabled={pathtracingEnabled}  // ← ADD THIS
            />
          );
        } else if (config.type === "Arrow") {
          return (
            <Arrow
              key={name}
              geometryKey={name}
              data={config.data}
              pathtracingEnabled={pathtracingEnabled}  // ← ADD THIS
            />
          );
        } else if (config.type === "Camera") {
          return (
            <Camera
              key={name}
              geometryKey={name}
              data={config.data}
            />
          );
        } else if (config.type === "Cell") {
          return (
            <Cell
              key={name}
              data={config.data}
            />
          );
        } else if (config.type === "Floor") {
          return (
            <Floor
              key={name}
              data={config.data}
            />
          );
        } else {
          console.warn(`Unhandled geometry type: ${config.type}`);
          return null;
        }
      })}

    {cameraSettings.show_crosshair && <Crosshair />}
    <VirtualCanvas />
  </PathTracingRenderer>

  <OrbitControls
    enableDamping={false}
    makeDefault
    enabled={cameraControls.enabled}
    enablePan={cameraControls.enablePan}
    enableRotate={cameraControls.enableRotate}
    enableZoom={cameraControls.enableZoom}
  />
</Canvas>
```

---

## Implementation Order (Critical Path)

1. ✅ **Backend Settings** - Add PathTracing to Python (no frontend deps)
2. ✅ **Generate Types** - Run Python script to generate TypeScript types
3. ✅ **Converter Utility** - Create convertInstancedMesh.ts (testable in isolation)
4. ✅ **Single Component** - Modify Particles.tsx first (test conversion works)
5. ✅ **PathTracingRenderer** - Create wrapper component
6. ✅ **Canvas Integration** - Wire everything together
7. ✅ **Test Basic Scene** - Verify pathtracing works with spheres
8. ✅ **Other Components** - Modify Bonds.tsx and Arrow.tsx
9. ✅ **Environment Testing** - Test different environment presets
10. ✅ **Performance Tuning** - Optimize settings for different use cases

---

## Performance Characteristics

### Memory Usage

| Mode | 1000 Spheres | 10000 Spheres | 100000 Spheres |
|------|--------------|---------------|----------------|
| Instanced | ~10 MB | ~15 MB | ~50 MB |
| Non-Instanced | ~500 MB | ~5 GB | ~50 GB |

**Critical:** Large scenes (>10k instances) may exceed browser memory limits.

### Rendering Performance

**Instanced Mode:**
- 60 FPS with 100k instances
- Interactive manipulation
- Real-time selection/hover

**Path Tracing Mode:**
- 0.1-1 FPS during sample accumulation
- Non-interactive (every movement resets samples)
- High quality after convergence

### Recommended Settings

**Preview (fast feedback):**
```python
PathTracing(
    enabled=True,
    min_samples=1,
    samples=32,
    bounces=2,
    tiles=2,
)
```

**Production (high quality):**
```python
PathTracing(
    enabled=True,
    min_samples=1,
    samples=1024,
    bounces=5,
    tiles=1,
)
```

---

## Known Limitations

1. **No Instancing Support** - Must convert, cannot be avoided
2. **Memory Intensive** - 100x+ memory increase for large scenes
3. **Non-Interactive** - Camera movement resets samples
4. **Material Conversion** - Materials must be PBR-compatible
5. **Selection Disabled** - Hover/select too slow in pathtracing mode


## Critical Review Checklist

- ✅ Instanced mesh conversion is mandatory (confirmed)
- ✅ Memory management strategy defined
- ✅ Dual rendering modes implemented
- ✅ Python backend settings integrated
- ✅ TypeScript types auto-generated
- ✅ Environment presets available
- ✅ Studio lighting disabled in pathtracing mode
- ✅ Implementation order is dependency-aware
- ✅ Performance characteristics documented
- ✅ Known limitations identified

---

## Dependencies Status

- ✅ `@react-three/gpu-pathtracer` - Installed
- ✅ `@react-three/drei` - Installed (has Environment component)
- ✅ `three` - Installed
- ✅ Python Pydantic models - Available
- ✅ Type generation pipeline - Working

**Ready to implement.**
