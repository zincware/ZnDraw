import { useFrame } from "@react-three/fiber";
import { usePathtracer } from "@react-three/gpu-pathtracer";
import { useEffect, useRef } from "react";
import { useAppStore } from "../store";
import type { PathTracing } from "../types/room-config";

/**
 * Component that watches for pathtracing update requests and calls
 * pathtracer.update() to rebuild the BVH acceleration structure.
 *
 * Uses useFrame to ensure update() is called after R3F has finished
 * processing all scene updates in the render loop.
 */
export function PathtracingUpdater({ settings }: { settings: PathTracing }) {
	const pathtracingNeedsUpdate = useAppStore(
		(state) => state.pathtracingNeedsUpdate,
	);
	const clearPathtracingUpdate = useAppStore(
		(state) => state.clearPathtracingUpdate,
	);
	const currentFrame = useAppStore((state) => state.currentFrame);

	const { update } = usePathtracer();

	// Track if this is the first render to avoid calling update on mount
	const isFirstRender = useRef(true);
	// Pending update flag - processed in useFrame to ensure scene is stable
	const pendingUpdate = useRef(false);

	// Capture manual update requests
	useEffect(() => {
		if (pathtracingNeedsUpdate) {
			pendingUpdate.current = true;
			clearPathtracingUpdate();
		}
	}, [pathtracingNeedsUpdate, clearPathtracingUpdate]);

	// Capture frame change updates
	useEffect(() => {
		if (isFirstRender.current) {
			isFirstRender.current = false;
			return;
		}
		pendingUpdate.current = true;
	}, [currentFrame]);

	// Capture settings change updates
	useEffect(() => {
		if (isFirstRender.current) {
			return;
		}
		pendingUpdate.current = true;
	}, [
		settings.samples,
		settings.min_samples,
		settings.bounces,
		settings.tiles,
		settings.environment_preset,
		settings.environment_intensity,
		settings.environment_blur,
		settings.environment_background,
	]);

	// Execute update in useFrame - guarantees R3F has processed scene updates
	useFrame(() => {
		if (pendingUpdate.current && update) {
			update();
			pendingUpdate.current = false;
		}
	});

	return null;
}
