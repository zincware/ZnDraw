import { Environment } from "@react-three/drei";
import { Pathtracer } from "@react-three/gpu-pathtracer";
import type { ReactNode } from "react";
import { PathtracingUpdater } from "./PathtracingUpdater";
import { PathtracingCaptureProvider } from "./three/PathtracingScreenshotCapture";

/**
 * PathTracing geometry data from backend.
 * Mirrors the Pydantic PathTracing model.
 */
interface PathTracingData {
	active?: boolean;
	min_samples?: number;
	samples?: number;
	bounces?: number;
	tiles?: number;
	environment_preset?: string;
	environment_intensity?: number;
	environment_blur?: number;
	environment_background?: boolean;
}

interface PathTracingRendererProps {
	pathtracingData?: PathTracingData;
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
	pathtracingData,
	children,
}: PathTracingRendererProps) {
	const {
		active = false,
		min_samples = 1,
		samples = 256,
		bounces = 3,
		tiles = 1,
		environment_preset = "studio",
		environment_intensity = 1.0,
		environment_blur = 0.0,
		environment_background = false,
	} = pathtracingData ?? {};

	return (
		<Pathtracer
			minSamples={min_samples}
			samples={samples}
			bounces={bounces}
			tiles={tiles}
			enabled={active}
		>
			{active && <PathtracingUpdater pathtracingData={pathtracingData!} />}
			{active && <PathtracingCaptureProvider />}

			{/* Environment lighting - only when pathtracing enabled */}
			{active && environment_preset !== "none" && (
				<Environment
					preset={
						environment_preset as
							| "apartment"
							| "city"
							| "dawn"
							| "forest"
							| "lobby"
							| "night"
							| "park"
							| "studio"
							| "sunset"
							| "warehouse"
					}
					background={environment_background}
					backgroundBlurriness={environment_blur}
					environmentIntensity={environment_intensity}
				/>
			)}

			{children}
		</Pathtracer>
	);
}
