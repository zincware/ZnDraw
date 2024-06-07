import { useState, useEffect, useRef, useMemo } from "react";
import * as THREE from "three";

import { useFrame, useThree } from "@react-three/fiber";
import { Line } from "@react-three/drei";

export interface Frame {
  arrays: { colors: Array<string>; radii: Array<number> };
  calc: any;
  cell: number[][];
  connectivity: Array<[number, number, number]>;
  info: any;
  numbers: number[];
  pbc: boolean[];
  positions: Array<[number, number, number]>;
}
``;

export interface Frames {
  [key: number]: Frame;
}

export const Player = ({
  playing,
  step,
  setStep,
  fps,
  length,
  loop,
  togglePlaying: setPlaying,
}: any) => {
  useFrame(({ clock }) => {
    const a = clock.getElapsedTime();
    if (a > 1 / fps) {
      if (playing) {
        if (step === length - 1) {
          if (!loop) {
            setPlaying(!playing);
          } else {
            setStep(0);
          }
        } else {
          setStep(step + 1);
        }
      }
      // reset the clock
      clock.start();
    }
  });
  return null;
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
}) => {
  const meshRef = useRef();

  const count = Object.keys(frame.positions).length;

  const { camera } = useThree();

  // TODO: look at the COM of the particles
  useEffect(() => {
    if (camera) {
      camera.position.set(-10, -10, -10); // Set initial position
      // camera.lookAt(0, 0, 0); // Make the camera look at (5, 5, 5)
    }
  }, [camera]);

  // attach key handler on press "c"
  useEffect(() => {
    const handleKeyDown = (event: KeyboardEvent) => {
      if (event.key === "c") {
        // get the center of mass of the selected particles
        const com = new THREE.Vector3();
        if (selectedIds.size === 0) {
          // get center from all particles
          frame.positions.forEach((position) => {
            com.add(new THREE.Vector3(...position));
          });
          com.divideScalar(frame.positions.length);
        } else {
          selectedIds.forEach((id) => {
            com.add(new THREE.Vector3(...frame.positions[id]));
          });
          com.divideScalar(selectedIds.size);
        }
        setOrbitControlsTarget([com.x, com.y, com.z]);
      }
    };
    window.addEventListener("keydown", handleKeyDown);
    return () => {
      window.removeEventListener("keydown", handleKeyDown);
    };
  }, [selectedIds, frame.positions, camera]);

  useEffect(() => {
    if (meshRef.current && count > 0) {
      const color = new THREE.Color();

      frame.positions.forEach((position, i) => {
        const matrix = new THREE.Matrix4().setPosition(...position);
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
          color.setHex(0xffa500);
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

  useEffect(() => {
    if (meshRef.current && count > 0) {
      const color = new THREE.Color();
      frame.positions.forEach((_, i) => {
        if (selectedIds.has(i)) {
          color.setHex(0xffa500);
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
          color.setHex(0xffa500);
          meshRef.current.setColorAt(i, color);
        } else {
          color.set(frame.arrays.colors[i]);
          meshRef.current.setColorAt(i, color);
        }
      });
      meshRef.current.instanceColor.needsUpdate = true;
    }
  });

  return (
    <instancedMesh
      ref={meshRef}
      args={[null, null, count]}
      onPointerOver={handlePointerOver}
      onPointerOut={handlePointerOut}
      onClick={handleClick}
      onPointerMove={handlePointerMove}
      castShadow
      receiveShadow
    >
      <sphereGeometry args={[1, 30, 30]} />
      {/* <meshPhongMaterial attach="material" shininess={100} /> */}
      <meshPhysicalMaterial
        attach="material"
        roughness={0.3}
        reflectivity={0.4}
      />
    </instancedMesh>
  );
};

export const BondInstances = ({
  frame,
  selectedIds,
  setSelectedIds,
}: {
  frame: Frame;
  selectedIds: Set<number>;
  setSelectedIds: any;
}) => {
  const meshRef = useRef<THREE.InstancedMesh | null>(null);

  const count = Object.keys(frame.connectivity).length * 2;

  const geometry = useMemo(() => {
    const _geometry = new THREE.CylinderGeometry(0.14, 0.14, 1, 30, 30);
    _geometry.translate(0, 0.25, 0);
    _geometry.rotateX(Math.PI / 2);
    _geometry.scale(1, 1, 0.5);
    return _geometry;
  }, []);

  useEffect(() => {
    if (meshRef.current && count > 0) {
      const color = new THREE.Color();
      const posA = new THREE.Vector3();
      const posB = new THREE.Vector3();
      const matrix = new THREE.Matrix4();
      const up = new THREE.Vector3(0, 1, 0);
      const direction = new THREE.Vector3();
      const scale = new THREE.Vector3(1, 1, 1);

      frame.connectivity.forEach(([a, b, order], i) => {
        posA.set(
          frame.positions[a][0],
          frame.positions[a][1],
          frame.positions[a][2],
        );
        posB.set(
          frame.positions[b][0],
          frame.positions[b][1],
          frame.positions[b][2],
        );

        direction.subVectors(posB, posA).normalize();
        const distance = posA.distanceTo(posB);
        scale.set(1, 1, distance);

        matrix.lookAt(posA, posB, up);
        matrix.scale(scale);
        matrix.setPosition(posA.lerp(posB, 0.5));
        meshRef.current.setMatrixAt(i * 2, matrix);
        if (selectedIds.has(a)) {
          color.setHex(0xffa500);
        } else {
          color.set(frame.arrays.colors[a]);
        }
        meshRef.current.setColorAt(i * 2, color);

        matrix.lookAt(posB, posA, up);
        matrix.scale(scale);
        matrix.setPosition(posB.lerp(posA, 0.5));
        meshRef.current.setMatrixAt(i * 2 + 1, matrix);
        if (selectedIds.has(b)) {
          color.setHex(0xffa500);
        } else {
          color.set(frame.arrays.colors[b]);
        }
        meshRef.current.setColorAt(i * 2 + 1, color);
      });

      meshRef.current.instanceMatrix.needsUpdate = true;
      meshRef.current.instanceColor.needsUpdate = true;
    }
  }, [frame]);

  // if selectedIds changes, update the color of the bonds corresponding to the selected particles
  useEffect(() => {
    if (meshRef.current && count > 0) {
      const color = new THREE.Color();
      frame.connectivity.forEach(([a, b], i) => {
        if (selectedIds.has(a)) {
          color.setHex(0xffa500);
          meshRef.current.setColorAt(i * 2, color);
        } else if (selectedIds.has(b)) {
          color.setHex(0xffa500);
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
  }, [selectedIds]);

  return (
    <instancedMesh
      ref={meshRef}
      args={[geometry, null, count]}
      castShadow
      receiveShadow
    >
      {/* <meshPhongMaterial attach="material" shininess={100} /> */}
      <meshPhysicalMaterial
        attach="material"
        roughness={0.3}
        reflectivity={0.4}
      />
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