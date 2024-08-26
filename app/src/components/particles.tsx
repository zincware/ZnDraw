import React, { useState, useEffect, useRef, useMemo } from "react";
import * as THREE from "three";

import { useFrame, useThree } from "@react-three/fiber";
import { Line } from "@react-three/drei";
import Arrows from "./meshes";
import { IndicesState } from "./utils";

export interface Frame {
  arrays: { colors: Array<string>; radii: Array<number> };
  calc: any;
  cell: number[][];
  connectivity: Array<[number, number, number]>;
  info: any;
  numbers: number[];
  pbc: boolean[];
  positions: THREE.Vector3[]; // only number[][] | before being mapped immediately
  vectors: [number, number, number][][];
}

export interface Frames {
  [key: number]: { _type: string; value: Frame };
}

interface PlayerProps {
  playing: boolean;
  step: number;
  setStep: (step: number) => void;
  fps: number;
  length: number;
  loop: boolean;
  togglePlaying: (playing: boolean) => void;
  selectedFrames: IndicesState;
}

export const Player = ({
  playing,
  step,
  setStep,
  fps,
  length,
  loop,
  togglePlaying: setPlaying,
  selectedFrames,
}: PlayerProps) => {
  useFrame(({ clock }) => {
    const a = clock.getElapsedTime();
    if (a > 1 / fps) {
      if (playing) {
        if (step === length - 1) {
          if (!loop) {
            setPlaying(!playing);
          } else {
            if ((selectedFrames.indices.size > 0) && selectedFrames.active) {
              setStep(Math.min(...selectedFrames.indices));
            } else {
              setStep(0);
            }
          }
        } else {
          if ((selectedFrames.indices.size > 0) && selectedFrames.active) {
            const nextFrame = Array.from(selectedFrames.indices).find(
              (frame) => frame > step,
            );
            if (nextFrame) {
              setStep(nextFrame);
            } else {
              setStep(Math.min(...selectedFrames.indices));
            }
          } else {
            setStep(step + 1);
          }
        }
      }
      // reset the clock
      clock.start();
    }
  });
  return null;
};

export const getCentroid = (
  positions: THREE.Vector3[],
  selection: Set<number>,
) => {
  if (selection.size > 0) {
    const centroid = new THREE.Vector3();
    selection.forEach((i) => {
      centroid.add(positions[i]);
    });
    centroid.divideScalar(selection.size);
    return centroid;
  } else {
    const centroid = new THREE.Vector3();
    positions.forEach((position) => {
      centroid.add(position);
    });
    centroid.divideScalar(positions.length);
    return centroid;
  }
};

export const ParticleInstances = ({
  frame,
  selectedIds,
  setSelectedIds,
  isDrawing,
  points,
  setPoints,
  setOrbitControlsTarget,
  hoveredId,
  setHoveredId,
  setTriggerSelection,
  sceneSettings,
}: {
  frame: Frame;
  selectedIds: Set<number>;
  setSelectedIds: any;
  isDrawing: boolean;
  points: THREE.Vector3[];
  setPoints: any;
  setOrbitControlsTarget: any;
  hoveredId: number | null;
  setHoveredId: any;
  setTriggerSelection: any;
  sceneSettings: any;
}) => {
  const meshRef = useRef();

  const count = Object.keys(frame.positions).length;
  const originalScale = useRef<number>(1);
  const sphereRef = useRef<THREE.InstancedMesh | null>(null);

  const { camera } = useThree();

  // TODO: look at the COM of the particles
  useEffect(() => {
    if (camera) {
      camera.position.set(-10, -10, -10); // Set initial position
      // camera.lookAt(0, 0, 0); // Make the camera look at (5, 5, 5)
    }
  }, [camera]);

  useEffect(() => {
    if (meshRef.current && count > 0) {
      const color = new THREE.Color();

      frame.positions.forEach((position, i) => {
        const matrix = new THREE.Matrix4().setPosition(position);
        matrix.scale(
          new THREE.Vector3(
            frame.arrays.radii[i],
            frame.arrays.radii[i],
            frame.arrays.radii[i],
          ),
        );
        meshRef.current.setMatrixAt(i, matrix);
        // update color
        if (selectedIds.has(i)) {
          color.setStyle(sceneSettings.selection_color);
        } else {
          color.set(frame.arrays.colors[i]);
        }
        meshRef.current.setColorAt(i, color);
      });
      meshRef.current.instanceMatrix.needsUpdate = true;
      meshRef.current.instanceColor.needsUpdate = true;
    }
  }, [frame.positions, frame.arrays.colors, frame.arrays.radii]);

  const handlePointerOver = (event) => {
    event.stopPropagation();
    setHoveredId(event.instanceId);
    // // detect shift and control key being pressed at the same time
    if (event.shiftKey && event.ctrlKey) {
      selectedIds.add(event.instanceId);
      setSelectedIds(new Set(selectedIds));
    }
  };

  const handlePointerMove = (event) => {
    event.stopPropagation();
    if (isDrawing) {
      // replace the last point with the hovered point
      let i = 0;
      while (
        i < event.intersections.length &&
        !event.intersections[i].object.visible
      ) {
        i++;
      }

      setPoints((prevPoints: THREE.Vector3[]) => [
        ...prevPoints.slice(0, prevPoints.length - 1),
        event.intersections[i].point,
      ]);
    }
  };

  const handlePointerOut = (event) => {
    event.stopPropagation();
    setHoveredId(null);
  };

  const handleClick = (event) => {
    // if not shift key, the selected particle is the only one selected, if already selected, no particle is selected
    // if shift key, the selected particle is added to the selected particles, if already selected, it is removed
    if (event.detail !== 1) return; // prevent double click tirgger
    if (!event.shiftKey) {
      if (selectedIds.has(event.instanceId)) {
        setSelectedIds(new Set());
      } else {
        setSelectedIds(new Set([event.instanceId]));
      }
    } else {
      if (selectedIds.has(event.instanceId)) {
        selectedIds.delete(event.instanceId);
        setSelectedIds(new Set(selectedIds));
      } else {
        selectedIds.add(event.instanceId);
        setSelectedIds(new Set(selectedIds));
      }
    }
    event.stopPropagation();
  };

  const handleDoubleClick = (event) => {
    setTriggerSelection(true);
    event.stopPropagation();
  };

  useEffect(() => {
    if (meshRef.current && count > 0) {
      const color = new THREE.Color();
      frame.positions.forEach((_, i) => {
        if (selectedIds.has(i)) {
          color.setStyle(sceneSettings.selection_color);
          meshRef.current.setColorAt(i, color);
        }
      });
      meshRef.current.instanceColor.needsUpdate = true;
    }
  }, [selectedIds]);

  // useFrame to change color if hovered
  useFrame(() => {
    if (meshRef.current && count > 0) {
      const color = new THREE.Color();
      frame.positions.forEach((_, i) => {
        if (i === hoveredId) {
          color.setHex(0xf05000);
          meshRef.current.setColorAt(i, color);
        } else if (selectedIds.has(i)) {
          color.setStyle(sceneSettings.selection_color);
          meshRef.current.setColorAt(i, color);
        } else {
          color.set(frame.arrays.colors[i]);
          meshRef.current.setColorAt(i, color);
        }
      });
      meshRef.current.instanceColor.needsUpdate = true;
    }
  });

  useEffect(() => {
    if (sphereRef.current) {
      sphereRef.current.scale(
        1 / originalScale.current,
        1 / originalScale.current,
        1 / originalScale.current,
      );
      originalScale.current = sceneSettings.particle_size;
      sphereRef.current.scale(
        sceneSettings.particle_size,
        sceneSettings.particle_size,
        sceneSettings.particle_size,
      );
    }
  }, [sceneSettings.particle_size]);

  return (
    <instancedMesh
      ref={meshRef}
      args={[null, null, count]}
      onPointerOver={handlePointerOver}
      onPointerOut={handlePointerOut}
      onClick={handleClick}
      onDoubleClick={handleDoubleClick}
      onPointerMove={handlePointerMove}
      castShadow
      receiveShadow
      frustumCulled={false}
    >
      <sphereGeometry args={[1, 30, 30]} ref={sphereRef} />
      {sceneSettings.material === "MeshPhysicalMaterial" && (
        <meshPhysicalMaterial
          attach="material"
          roughness={0.3}
          reflectivity={0.4}
        />
      )}
      {sceneSettings.material === "MeshToonMaterial" && (
        <meshToonMaterial attach="material" />
      )}
      {sceneSettings.material === "MeshStandardMaterial" && (
        <meshStandardMaterial attach="material" />
      )}
      {sceneSettings.material === "MeshBasicMaterial" && (
        <meshBasicMaterial attach="material" />
      )}
    </instancedMesh>
  );
};

export const BondInstances = ({
  frame,
  selectedIds,
  hoveredId,
  sceneSettings,
}: {
  frame: Frame;
  selectedIds: Set<number>;
  hoveredId: number;
  sceneSettings: any;
}) => {
  const meshRef = useRef<THREE.InstancedMesh | null>(null);

  const count = Object.keys(frame.connectivity).length * 2;
  const originalScale = useRef<number>(1);

  const geometry = useMemo(() => {
    const _geometry = new THREE.CylinderGeometry(0.14, 0.14, 1, 30, 30);
    _geometry.translate(0, 0.5, 0);
    _geometry.rotateX(Math.PI / 2);
    _geometry.scale(1, 1, 0.5);
    return _geometry;
  }, []);

  useEffect(() => {
    if (meshRef.current && count > 0) {
      const color = new THREE.Color();
      const matrix = new THREE.Matrix4();
      const up = new THREE.Vector3(0, 1, 0);
      const direction = new THREE.Vector3();
      const scale = new THREE.Vector3(1, 1, 1);

      const createTransformationMatrix = (
        posA: THREE.Vector3,
        posB: THREE.Vector3,
      ) => {
        direction.subVectors(posB, posA).normalize();
        const distance = posA.distanceTo(posB);
        scale.set(1, 1, distance);

        matrix.lookAt(posA, posB, up);
        matrix.scale(scale);
        matrix.setPosition(posA.clone().lerp(posB, 0.5));

        return matrix.clone(); // Clone to avoid overwriting
      };

      frame.connectivity.forEach(([a, b], i) => {
        const posA = frame.positions[a] as THREE.Vector3;
        const posB = frame.positions[b] as THREE.Vector3;
        if (!posA || !posB) {
          console.error("Connected particles not found");
          return;
        }

        // Set matrix and color for the bond from A to B
        meshRef.current!.setMatrixAt(
          i * 2,
          createTransformationMatrix(posA, posB),
        );
        color.setStyle(
          selectedIds.has(a)
            ? sceneSettings.selection_color
            : frame.arrays.colors[a],
        );
        meshRef.current!.setColorAt(i * 2, color);

        // Set matrix and color for the bond from B to A
        meshRef.current!.setMatrixAt(
          i * 2 + 1,
          createTransformationMatrix(posB, posA),
        );
        color.setStyle(
          selectedIds.has(b)
            ? sceneSettings.selection_color
            : frame.arrays.colors[b],
        );
        meshRef.current!.setColorAt(i * 2 + 1, color);
      });

      meshRef.current.instanceMatrix.needsUpdate = true;
      meshRef.current.instanceColor.needsUpdate = true;
    }
  }, [frame, selectedIds, count]);

  // if selectedIds changes, update the color of the bonds corresponding to the selected particles
  useEffect(() => {
    if (meshRef.current && count > 0) {
      const color = new THREE.Color();
      frame.connectivity.forEach(([a, b], i) => {
        if (selectedIds.has(a)) {
          color.setStyle(sceneSettings.selection_color);
          meshRef.current.setColorAt(i * 2, color);
        } else if (selectedIds.has(b)) {
          // TODO: check if this is required?
          color.setStyle(sceneSettings.selection_color);
          meshRef.current.setColorAt(i * 2 + 1, color);
        } else {
          color.set(frame.arrays.colors[a]);
          meshRef.current.setColorAt(i * 2, color);
          color.set(frame.arrays.colors[b]);
          meshRef.current.setColorAt(i * 2 + 1, color);
        }
      });
      meshRef.current.instanceColor.needsUpdate = true;
    }
  }, [selectedIds, sceneSettings]);

  useEffect(() => {
    if (meshRef.current) {
      geometry.scale(1 / originalScale.current, 1 / originalScale.current, 1);
      originalScale.current = sceneSettings.bond_size;
      geometry.scale(sceneSettings.bond_size, sceneSettings.bond_size, 1);
    }
  }, [sceneSettings.bond_size]);

  return (
    <instancedMesh
      ref={meshRef}
      args={[geometry, null, count]}
      castShadow
      receiveShadow
    >
      {sceneSettings.material === "MeshPhysicalMaterial" && (
        <meshPhysicalMaterial
          attach="material"
          roughness={0.3}
          reflectivity={0.4}
        />
      )}
      {sceneSettings.material === "MeshToonMaterial" && (
        <meshToonMaterial attach="material" />
      )}
      {sceneSettings.material === "MeshStandardMaterial" && (
        <meshStandardMaterial attach="material" />
      )}
      {sceneSettings.material === "MeshBasicMaterial" && (
        <meshBasicMaterial attach="material" />
      )}
    </instancedMesh>
  );
};

export const SimulationCell = ({
  frame,
  colorMode,
}: {
  frame: Frame;
  colorMode: string;
}) => {
  const [lineColor, setLineColor] = useState("black");

  useEffect(() => {
    if (colorMode === "light") {
      setLineColor("black");
    } else {
      setLineColor("white");
    }
  }, [colorMode]);

  const vertices = useMemo(() => {
    const cell = frame.cell;
    if (cell.length !== 3) {
      console.error("Invalid cell dimensions");
      return;
    }

    const origin = new THREE.Vector3(0, 0, 0);
    const vectors = cell.map((row) => new THREE.Vector3(...row));

    // Create the vertices of the box
    const v = [
      origin,
      vectors[0],
      vectors[1],
      vectors[0].clone().add(vectors[1]),
      vectors[2],
      vectors[0].clone().add(vectors[2]),
      vectors[1].clone().add(vectors[2]),
      vectors[0].clone().add(vectors[1]).add(vectors[2]),
    ];

    return [
      [
        // Bottom face
        v[0],
        v[1],
        v[1],
        v[3],
        v[3],
        v[2],
        v[2],
        v[0],
      ],
      [
        // Top face
        v[4],
        v[5],
        v[5],
        v[7],
        v[7],
        v[6],
        v[6],
        v[4],
      ],
      // Vertical lines

      [v[0], v[4]],
      [v[1], v[5]],
      [v[2], v[6]],
      [v[3], v[7]],
    ];
  }, [frame.cell]);

  return (
    <>
      {vertices && (
        <>
          <Line points={vertices[0]} color={lineColor} lineWidth={2} />
          <Line points={vertices[1]} color={lineColor} lineWidth={2} />
          <Line points={vertices[2]} color={lineColor} lineWidth={2} />
          <Line points={vertices[3]} color={lineColor} lineWidth={2} />
          <Line points={vertices[4]} color={lineColor} lineWidth={2} />
          <Line points={vertices[5]} color={lineColor} lineWidth={2} />
        </>
      )}
    </>
  );
};

interface PerParticleVectorsProps {
  frame: Frame | undefined;
  property: string;
  colorMode: string;
  arrowsConfig: any;
}

export const PerParticleVectors: React.FC<PerParticleVectorsProps> = ({
  frame,
  property,
  colorMode,
  arrowsConfig,
}) => {
  const [vectors, setVectors] = useState<
    { start: THREE.Vector3; end: THREE.Vector3 }[]
  >([]);

  const LineColor = colorMode === "light" ? "#454b66" : "#f5fdc6";
  const LineWidth = 2;
  const LineScale = 1;

  const [colorRange, setColorRange] = useState<[number, number]>(
    arrowsConfig.colorrange,
  );

  useEffect(() => {
    if (arrowsConfig.normalize) {
      const max = Math.max(
        ...vectors.map((vector) => vector.start.distanceTo(vector.end)),
      );
      setColorRange([0, max]);
    } else {
      setColorRange(arrowsConfig.colorrange);
    }
  }, [vectors, arrowsConfig.normalize, arrowsConfig.colorrange]);

  useEffect(() => {
    if (!frame || !frame.calc || !frame.calc[property]) {
      console.log(`Property ${property} not found in frame.calc`);
      setVectors([]);
      return;
    } else {
      console.log(`Property ${property} found in frame.calc`);
      const calculatedVectors = frame.calc[property].map((vector, i) => {
        const start = frame.positions[i];
        const end = start
          .clone()
          .add(
            new THREE.Vector3(vector[0], vector[1], vector[2]).multiplyScalar(
              LineScale,
            ),
          );
        return { start, end };
      });
      setVectors(calculatedVectors);
    }
  }, [frame, property]);

  return (
    <>
      <Arrows
        start={vectors.map((vec) => vec.start.toArray())}
        end={vectors.map((vec) => vec.end.toArray())}
        scale_vector_thickness={arrowsConfig.scale_vector_thickness}
        colormap={arrowsConfig.colormap}
        colorrange={colorRange}
        opacity={arrowsConfig.opacity}
        rescale={arrowsConfig.rescale}
      />
    </>
  );
};
