import { Pathtracer } from "@react-three/gpu-pathtracer";
import { Environment } from '@react-three/drei';
import type { PathTracing } from '../types/room-config';
import { PathtracingUpdater } from './PathtracingUpdater';

interface PathTracingRendererProps {
  settings?: PathTracing;
  children: React.ReactNode;
}

/**
 * Wraps scene with GPU path tracing renderer and environment lighting.
 * When disabled, passes children through without modification.
 *
 * IMPORTANT: When path tracing is enabled:
 * - All user interactions (click, hover, selection) are disabled
 * - Instanced meshes are converted to merged meshes (required for path tracer)
 * - Studio lighting is automatically disabled (environment provides lighting)
 * - Scene updates trigger BVH rebuild (can be expensive for large scenes)
 */
export function PathTracingRenderer({
  settings,
  children
}: PathTracingRendererProps) {
  // Handle undefined settings (loading state)
  if (!settings) {
    return <>{children}</>;
  }

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
      enabled={true}
    >
      {/* Pathtracing updater - watches for scene changes and calls update() */}
      <PathtracingUpdater settings={settings} />

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
