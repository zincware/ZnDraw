import { useState, useEffect } from "react";
import * as THREE from "three";
import type { HSLColor } from "../components/utils";
import { useFrameFetching } from "../components/api";
import type { Frame } from "../components/particles";

// Colormap definitions
const VIRIDIS_COLORMAP: HSLColor[] = [
	[0.65, 0.9, 0.15],  // Dark blue-purple
	[0.55, 0.85, 0.35],  // Blue-cyan
	[0.45, 0.9, 0.55],   // Green-cyan
	[0.25, 0.95, 0.75],  // Yellow-green
	[0.15, 1.0, 0.85],   // Yellow
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
	vectors: string[];
	colormap: string;
	vector_scale: number;
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
	perParticleVectors: { start: THREE.Vector3; end: THREE.Vector3 }[];
	vectorColormap: HSLColor[];
}

export const useVectorManager = ({
	token,
	step,
	currentFrame,
	vectorConfig,
}: UseVectorManagerProps): UseVectorManagerReturn => {
	const [vectorProperties, setVectorProperties] = useState<{ [key: string]: any }>({});
	const [perParticleVectors, setPerParticleVectors] = useState<{ start: THREE.Vector3; end: THREE.Vector3 }[]>([]);
	
	// Computed colormap based on vectorConfig
	const vectorColormap = getColormap(vectorConfig.colormap);

	// Custom hook for vector setup
	const { getFrameFromCon } = useFrameFetching(token);

	// Setup vectors when needed
	useEffect(() => {
		const setupVectorsFromFrame = async (frame: any, requestedProperties: string[]) => {
			if (!frame) {
				setVectorProperties({});
				return;
			}

			try {
				const [calc, arrays] = await Promise.all([frame.calc, frame.arrays]);

				const [calcKeys, arraysKeys] = await Promise.all([
					calc.keys(),
					arrays.keys(),
				]);

				const vectorPropertiesData: { [key: string]: any } = {};
				const allKeys = new Set([...calcKeys, ...arraysKeys]);

				// Only process requested properties
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
						console.warn(`Property '${property}' is not a valid 3D vector array`);
					}
				}

				setVectorProperties(vectorPropertiesData);
			} catch (error) {
				console.error("Error setting up vectors:", error);
				setVectorProperties({});
			}
		};

		if (token && vectorConfig.vectors && vectorConfig.vectors.length > 0) {
			console.log("Setting up vectors for properties:", vectorConfig.vectors);
			
			const loadVectors = async () => {
				try {
					const frame = await getFrameFromCon(step);
					await setupVectorsFromFrame(frame, vectorConfig.vectors);
				} catch (error) {
					console.error("Error loading vectors:", error);
					setVectorProperties({});
				}
			};

			loadVectors();
		} else {
			setVectorProperties({});
		}
	}, [token, step, vectorConfig.vectors, getFrameFromCon]);

	// Calculate per-particle vectors when vectorProperties or currentFrame changes
	useEffect(() => {
		if (!currentFrame || !vectorConfig.vectors || vectorConfig.vectors.length === 0) {
			setPerParticleVectors([]);
			return;
		}

		// Combine all selected vector properties
		const allVectors: { start: THREE.Vector3; end: THREE.Vector3 }[] = [];
		
		for (const vectorProperty of vectorConfig.vectors) {
			const frameData = vectorProperties[vectorProperty];
			
			if (frameData && frameData.length === currentFrame.positions.length) {
				const calculatedVectors = frameData.map(
					(vector: [number, number, number], i: number) => ({
						start: currentFrame.positions[i],
						end: currentFrame.positions[i].clone().add(
							new THREE.Vector3(...vector).multiplyScalar(vectorConfig.vector_scale)
						),
					})
				);
				allVectors.push(...calculatedVectors);
			}
		}
		
		console.log("Calculated per-particle vectors:", allVectors);
		setPerParticleVectors(allVectors);
	}, [currentFrame, vectorProperties, vectorConfig.vectors, vectorConfig.vector_scale]);

	return {
		vectorProperties,
		setVectorProperties,
		perParticleVectors,
		vectorColormap,
	};
};