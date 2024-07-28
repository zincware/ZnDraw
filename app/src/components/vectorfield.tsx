import React from "react";
import Arrow from "./meshes";

interface VectorFieldProps {
  vectors: [number, number, number][][];
  showArrows?: boolean;
}

const VectorField: React.FC<VectorFieldProps> = ({
  vectors,
  showArrows = true,
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
            />
          )}
        </React.Fragment>
      ))}
    </>
  );
};

export default VectorField;
