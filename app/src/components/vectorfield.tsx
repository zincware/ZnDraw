import type React from "react";
import { useEffect, useMemo, useState } from "react";
import * as THREE from "three";
import Arrows from "./meshes";
import type { HSLColor } from "./utils";

interface VectorFieldProps {
	vectors: [number, number, number][][];
	showArrows?: boolean;
	pathTracingSettings: any | undefined;
	arrowsConfig: {
		normalize: boolean;
		scale_vector_thickness: boolean;
		colormap: HSLColor[];
		colorrange: [number, number];
		opacity: number;
	};
}

export const VectorField: React.FC<VectorFieldProps> = ({
	vectors,
	showArrows = true,
	arrowsConfig,
	pathTracingSettings,
}) => {
	const [colorRange, setColorRange] = useState<[number, number]>(
		Array.isArray(arrowsConfig.colorrange) ? arrowsConfig.colorrange : [0, 1],
	);

	useEffect(() => {
		if (vectors.length === 0) {
			setColorRange([0, 1]);
			return;
		}

		if (arrowsConfig.normalize) {
			const max = Math.max(
				...vectors.map((vector) =>
					new THREE.Vector3(...vector[0]).distanceTo(
						new THREE.Vector3(...vector[1]),
					),
				),
			);
			setColorRange([0, max]);
		} else {
			setColorRange(
				Array.isArray(arrowsConfig.colorrange)
					? arrowsConfig.colorrange
					: [0, 1],
			);
		}
	}, [vectors, arrowsConfig.normalize, arrowsConfig.colorrange]);

	const startMap = useMemo(() => vectors.map((vector) => vector[0]), [vectors]);
	const endMap = useMemo(() => vectors.map((vector) => vector[1]), [vectors]);

	if (vectors.length === 0) {
		return null;
	}

	return (
		<Arrows
			start={startMap}
			end={endMap}
			scale_vector_thickness={arrowsConfig.scale_vector_thickness}
			colormap={arrowsConfig.colormap}
			colorrange={colorRange}
			opacity={arrowsConfig.opacity}
			rescale={arrowsConfig.rescale}
			pathTracingSettings={pathTracingSettings}
		/>
	);
};

export default VectorField;
