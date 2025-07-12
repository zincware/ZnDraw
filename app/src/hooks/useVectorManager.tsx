import { useState, useEffect } from "react";
import * as THREE from "three";
import type { HSLColor } from "../components/utils";
import { useFrameFetching } from "../components/api";
import type { Frame } from "../components/particles";

// Colormap definitions
const VIRIDIS_COLORMAP: HSLColor[] = [
	[0.65, 0.9, 0.15], // Dark blue-purple
	[0.55, 0.85, 0.35], // Blue-cyan
	[0.45, 0.9, 0.55], // Green-cyan
	[0.25, 0.95, 0.75], // Yellow-green
	[0.15, 1.0, 0.85], // Yellow
];

const getColormap = (colormapName: string): HSLColor[] => {
	switch (colormapName) {
		case "viridis":
			return VIRIDIS_COLORMAP;
		default:
			return VIRIDIS_COLORMAP; // Default to viridis
	}
};

interface VectorConfig {
	vectors?: string[];
	vector_scale?: number;
	vectorfield?: boolean;
	vector_colors?: { [key: string]: string };
	[key: string]: any;
}

interface UseVectorManagerProps {
	token: string;
	step: number;
	currentFrame: Frame;
	vectorConfig: VectorConfig;
}

interface UseVectorManagerReturn {
	vectorProperties: { [key: string]: any };
	setVectorProperties: (properties: { [key: string]: any }) => void;
	perParticleVectors: {
		start: THREE.Vector3;
		end: THREE.Vector3;
		vectorType: string;
	}[];
	vectorFieldData: [number, number, number][][];
	vectorColormap: HSLColor[];
}

export const useVectorManager = ({
	token,
	step,
	currentFrame,
	vectorConfig,
}: UseVectorManagerProps): UseVectorManagerReturn => {
	const [vectorProperties, setVectorProperties] = useState<{
		[key: string]: any;
	}>({});
	const [perParticleVectors, setPerParticleVectors] = useState<
		{ start: THREE.Vector3; end: THREE.Vector3; vectorType: string }[]
	>([]);
	const [vectorFieldData, setVectorFieldData] = useState<
		[number, number, number][][]
	>([]);

	// Computed colormap based on vectorConfig (using default since colormap is no longer used)
	const vectorColormap = getColormap("viridis");

	// Custom hook for vector setup
	const { getFrameFromCon } = useFrameFetching(token);

	// Setup vectors when needed
	useEffect(() => {
		const setupVectorsFromFrame = async (
			frame: any,
			requestedProperties: string[],
		) => {
			if (!frame) {
				setVectorProperties({});
				setVectorFieldData([]);
				return;
			}

			try {
				const [calc, arrays, vectors] = await Promise.all([
					frame.calc,
					frame.arrays,
					frame.vectors,
				]);

				const [calcKeys, arraysKeys] = await Promise.all([
					calc.keys(),
					arrays.keys(),
				]);

				const vectorPropertiesData: { [key: string]: any } = {};
				const allKeys = new Set([...calcKeys, ...arraysKeys]);

				// Only process requested properties for per-particle vectors
				for (const property of requestedProperties) {
					if (!allKeys.has(property)) {
						console.warn(`Requested property '${property}' not found in frame`);
						continue;
					}

					let data;
					if (calcKeys.includes(property)) {
						data = await calc[property];
					} else if (arraysKeys.includes(property)) {
						data = await arrays[property];
					}

					if (
						data &&
						Array.isArray(data) &&
						data.length > 0 &&
						Array.isArray(data[0]) &&
						data[0].length === 3
					) {
						vectorPropertiesData[property] = data;
					} else {
						console.warn(
							`Property '${property}' is not a valid 3D vector array`,
						);
					}
				}

				setVectorProperties(vectorPropertiesData);

				// Load vector field data from info.vectors if vectorfield is enabled
				if (vectorConfig.vectorfield && vectors) {
					try {
						if (vectors && Array.isArray(vectors)) {
							setVectorFieldData(vectors);
						} else {
							setVectorFieldData([]);
						}
					} catch (error) {
						console.warn("No vector field data found in frame.info.vectors");
						setVectorFieldData([]);
					}
				} else {
					setVectorFieldData([]);
				}
			} catch (error) {
				console.error("Error setting up vectors:", error);
				setVectorProperties({});
				setVectorFieldData([]);
			}
		};

		const shouldLoadVectors =
			token && vectorConfig.vectors && vectorConfig.vectors.length > 0;
		const shouldLoadVectorField = token && (vectorConfig.vectorfield || false);

		if (shouldLoadVectors || shouldLoadVectorField) {
			const loadVectors = async () => {
				try {
					const frame = await getFrameFromCon(step);
					await setupVectorsFromFrame(frame, vectorConfig.vectors || []);
				} catch (error) {
					console.error("Error loading vectors:", error);
					setVectorProperties({});
					setVectorFieldData([]);
				}
			};

			loadVectors();
		} else {
			setVectorProperties({});
			setVectorFieldData([]);
		}
	}, [
		token,
		step,
		vectorConfig.vectors,
		vectorConfig.vectorfield,
		getFrameFromCon,
	]);

	// Calculate per-particle vectors when vectorProperties or currentFrame changes
	useEffect(() => {
		if (
			!currentFrame ||
			!vectorConfig.vectors ||
			vectorConfig.vectors.length === 0
		) {
			setPerParticleVectors([]);
			return;
		}

		// Combine all selected vector properties with type information
		const allVectors: {
			start: THREE.Vector3;
			end: THREE.Vector3;
			vectorType: string;
		}[] = [];

		for (const vectorProperty of vectorConfig.vectors) {
			const frameData = vectorProperties[vectorProperty];

			if (frameData && frameData.length === currentFrame.positions.length) {
				const calculatedVectors = frameData.map(
					(vector: [number, number, number], i: number) => ({
						start: currentFrame.positions[i],
						end: currentFrame.positions[i]
							.clone()
							.add(
								new THREE.Vector3(...vector).multiplyScalar(
									vectorConfig.vector_scale || 1.0,
								),
							),
						vectorType: vectorProperty, // Add vector type information
					}),
				);
				allVectors.push(...calculatedVectors);
			}
		}

		console.log("Calculated per-particle vectors:", allVectors);
		setPerParticleVectors(allVectors);
	}, [
		currentFrame,
		vectorProperties,
		vectorConfig.vectors,
		vectorConfig.vector_scale,
	]);

	return {
		vectorProperties,
		setVectorProperties,
		perParticleVectors,
		vectorFieldData,
		vectorColormap,
	};
};
