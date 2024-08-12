import React, { useEffect, useState } from "react";
import * as THREE from "three";
import Arrows from "./meshes";
import { HSLColor } from "./utils";

interface VectorFieldProps {
  vectors: [number, number, number][][];
  showArrows?: boolean;
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
}) => {
  const [colorRange, setColorRange] = useState<[number, number]>(
    arrowsConfig.colorrange,
  );

  useEffect(() => {
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
      setColorRange(arrowsConfig.colorrange);
    }
  }, [vectors, arrowsConfig.normalize, arrowsConfig.colorrange]);

  return (
    <Arrows
      start={vectors.map((vector) => vector[0])}
      end={vectors.map((vector) => vector[1])}
      scale_vector_thickness={arrowsConfig.scale_vector_thickness}
      colormap={arrowsConfig.colormap}
      colorrange={colorRange}
      opacity={arrowsConfig.opacity}
      rescale={arrowsConfig.rescale}
    />
  );
};

export default VectorField;
