import React from "react";
import Arrow from "./meshes";

interface VectorFieldProps {
  vectors: [number, number, number][][];
  showArrows?: boolean;
  scale_vector_thickness?: boolean;
}

const VectorField: React.FC<VectorFieldProps> = ({
  vectors,
  showArrows = true,
  scale_vector_thickness = false,
}) => {
  return (
    <>
      {vectors.map((vector, index) => (
        <React.Fragment key={index}>
          {/* <Line points={vector} color="blue" lineWidth={2} /> */}
          {showArrows && vector.length > 1 && (
            <Arrow
              start={vector[vector.length - 2]}
              end={vector[vector.length - 1]}
              scale_vector_thickness={scale_vector_thickness}
            />
          )}
        </React.Fragment>
      ))}
    </>
  );
};

export default VectorField;
