import type { ReactNode } from "react";
import { Pathtracer } from "@react-three/gpu-pathtracer";
import { Environment } from "@react-three/drei";
import type { PathTracing } from "../types/room-config";
import { PathtracingUpdater } from "./PathtracingUpdater";
import { PathtracingCaptureProvider } from "./three/PathtracingScreenshotCapture";

interface PathTracingRendererProps {
	settings?: PathTracing;
	children: ReactNode;
}

/**
 * Wraps scene with GPU path tracing renderer and environment lighting.
 *
 * When path tracing is enabled:
 * - All user interactions (click, hover, selection) are disabled
 * - Instanced meshes are converted to merged meshes (required for path tracer)
 * - Studio lighting is automatically disabled (environment provides lighting)
 * - Scene updates trigger BVH rebuild (can be expensive for large scenes)
 */
export function PathTracingRenderer({
	settings,
	children,
}: PathTracingRendererProps) {
	const {
		enabled = false,
		min_samples = 1,
		samples = 256,
		bounces = 3,
		tiles = 1,
		environment_preset = "studio",
		environment_intensity = 1.0,
		environment_blur = 0.0,
		environment_background = false,
	} = settings ?? {};
	return (
		<Pathtracer
			minSamples={min_samples}
			samples={samples}
			bounces={bounces}
			tiles={tiles}
			enabled={enabled}
		>
			{enabled && <PathtracingUpdater settings={settings!} />}
			{enabled && <PathtracingCaptureProvider />}

			{/* Environment lighting - only when pathtracing enabled */}
			{enabled && environment_preset !== "none" && (
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
