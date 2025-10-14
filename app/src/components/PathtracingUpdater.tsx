import { useEffect, useRef } from 'react';
import { usePathtracer } from '@react-three/gpu-pathtracer';
import { useAppStore } from '../store';
import type { PathTracing } from '../types/room-config';

/**
 * Component that watches for pathtracing update requests and calls
 * pathtracer.update() to rebuild the BVH acceleration structure.
 *
 * This is needed when:
 * - Scene changes dynamically (e.g., particles switching from instanced to individual meshes)
 * - Current frame changes (particles move, arrows change, etc.)
 * - Pathtracing settings change
 */
export function PathtracingUpdater({ settings }: { settings: PathTracing }) {
  const pathtracingNeedsUpdate = useAppStore((state) => state.pathtracingNeedsUpdate);
  const clearPathtracingUpdate = useAppStore((state) => state.clearPathtracingUpdate);
  const currentFrame = useAppStore((state) => state.currentFrame);

  const { update } = usePathtracer();

  // Track if this is the first render to avoid calling update on mount
  const isFirstRender = useRef(true);

  // Handle manual update requests (e.g., when switching to individual meshes)
  useEffect(() => {
    if (pathtracingNeedsUpdate && update) {
      update();
      clearPathtracingUpdate();
    }
  }, [pathtracingNeedsUpdate, update, clearPathtracingUpdate]);

  // Handle automatic updates when frame changes
  useEffect(() => {
    // Skip the first render
    if (isFirstRender.current) {
      isFirstRender.current = false;
      return;
    }

    if (update) {
      update();
    }
  }, [currentFrame, update]);

  // Handle automatic updates when pathtracing settings change
  useEffect(() => {
    // Skip the first render
    if (isFirstRender.current) {
      return;
    }

    if (update) {
      update();
    }
  }, [
    settings.samples,
    settings.min_samples,
    settings.bounces,
    settings.tiles,
    settings.environment_preset,
    settings.environment_intensity,
    settings.environment_blur,
    settings.environment_background,
    update,
  ]);

  // This component doesn't render anything
  return null;
}
