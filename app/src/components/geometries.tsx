import {
  Box,
  Circle,
  Plane,
  Cone,
  Cylinder,
  Ring,
  Torus,
  TorusKnot,
  Sphere,
  Dodecahedron,
  Icosahedron,
  Octahedron,
  Tetrahedron,
  Outlines,
} from "@react-three/drei";

import * as THREE from "three";

import { useRef, useEffect } from "react";

export interface BaseGeometry {
  material: {
    color: string;
    opacity: number;
    wireframe: boolean;
    outlines: boolean;
  };
  position: [number, number, number];
  rotation: [number, number, number];
  scale: [number, number, number];
  discriminator: string;
}

export interface BoxGeometry extends BaseGeometry {
  discriminator: "Box";
  width: number;
  height: number;
  depth: number;
}

export interface CircleGeometry extends BaseGeometry {
  discriminator: "Circle";
  radius: number;
}

export interface PlaneGeometry extends BaseGeometry {
  discriminator: "Plane";
  width: number;
  height: number;
}

export interface ConeGeometry extends BaseGeometry {
  discriminator: "Cone";
  radius: number;
  height: number;
}

export interface CylinderGeometry extends BaseGeometry {
  discriminator: "Cylinder";
  radius_top: number;
  radius_bottom: number;
  height: number;
}

export interface RingGeometry extends BaseGeometry {
  discriminator: "Ring";
  inner_radius: number;
  outer_radius: number;
}

export interface TorusGeometry extends BaseGeometry {
  discriminator: "Torus";
  radius: number;
  tube: number;
}

export interface SphereGeometry extends BaseGeometry {
  discriminator: "Sphere";
  radius: number;
}

export interface DodecahedronGeometry extends BaseGeometry {
  discriminator: "Dodecahedron";
  radius: number;
}

export interface IcosahedronGeometry extends BaseGeometry {
  discriminator: "Icosahedron";
  radius: number;
}

export interface OctahedronGeometry extends BaseGeometry {
  discriminator: "Octahedron";
  radius: number;
}

export interface TetrahedronGeometry extends BaseGeometry {
  discriminator: "Tetrahedron";
  radius: number;
}

export interface TorusKnotGeometry extends BaseGeometry {
  discriminator: "TorusKnot";
  radius: number;
  tube: number;
}

export type Geometry =
  | BoxGeometry
  | CircleGeometry
  | PlaneGeometry
  | ConeGeometry
  | CylinderGeometry
  | RingGeometry
  | SphereGeometry
  | DodecahedronGeometry
  | IcosahedronGeometry
  | TetrahedronGeometry
  | TorusKnotGeometry
  | OctahedronGeometry
  | TorusGeometry;

function GeometryComponent({
  geometry,
  onPointerMove,
  onPointerOver,
  onPointerOut,
}: {
  geometry: Geometry;
  onPointerMove: any;
  onPointerOver: any;
  onPointerOut: any;
}) {
  switch (geometry.discriminator) {
    case "TorusKnot":
      return (
        <TorusKnot
          args={[geometry.radius, geometry.tube]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </TorusKnot>
      );
    case "Tetrahedron":
      return (
        <Tetrahedron
          args={[geometry.radius]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Tetrahedron>
      );
    case "Octahedron":
      return (
        <Octahedron
          args={[geometry.radius]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Octahedron>
      );
    case "Icosahedron":
      return (
        <Icosahedron
          args={[geometry.radius]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Icosahedron>
      );
    case "Dodecahedron":
      return (
        <Dodecahedron
          args={[geometry.radius]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Dodecahedron>
      );
    case "Sphere":
      return (
        <Sphere
          args={[geometry.radius, 32, 32]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Sphere>
      );

    case "Box":
      return (
        <Box
          args={[geometry.width, geometry.height, geometry.depth]}
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Box>
      );
    case "Circle":
      return (
        <Circle
          args={[geometry.radius, 32]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Circle>
      );
    case "Plane":
      return (
        <Plane
          args={[geometry.width, geometry.height]}
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Plane>
      );
    case "Cone":
      return (
        <Cone
          args={[geometry.radius, geometry.height, 32]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Cone>
      );
    case "Cylinder":
      return (
        <Cylinder
          args={[
            geometry.radius_top,
            geometry.radius_bottom,
            geometry.height,
            32,
          ]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Cylinder>
      );
    case "Ring":
      return (
        <Ring
          args={[geometry.inner_radius, geometry.outer_radius, 32]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Ring>
      );
    case "Torus":
      return (
        <Torus
          args={[geometry.radius, geometry.tube, 16, 100]} // Adjust the segments as needed
          position={geometry.position}
          rotation={geometry.rotation}
          scale={geometry.scale}
          onPointerMove={onPointerMove}
          onPointerOver={onPointerOver}
          onPointerOut={onPointerOut}
        >
          <meshStandardMaterial
            attach="material"
            color={geometry.material.color}
            opacity={geometry.material.opacity}
            wireframe={geometry.material.wireframe}
            transparent={geometry.material.opacity < 1.0}
          />
          {geometry.material.outlines && (
            <Outlines
              thickness={0.05}
              color={geometry.material.color}
              opacity={geometry.material.opacity}
            />
          )}
        </Torus>
      );
    default:
      return null;
  }
}

export function Geometries({
  geometries,
  isDrawing,
  setPoints,
  setHoveredId,
}: {
  geometries: Geometry[];
  isDrawing: boolean;
  setPoints: any;
  setHoveredId: any;
}) {
  const handlePointerMove = (event: any) => {
    event.stopPropagation();

    if (isDrawing) {
      // replace the last point with the hovered point
      // find the first object in event.intersections that is visible
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

  const handlePointerOver = (event: any) => {
    setHoveredId(event.object.uuid);
  };

  const handlePointerOut = (event: any) => {
    setHoveredId(null);
  };

  return (
    <>
      {geometries.map((geometry, index) => (
        <GeometryComponent
          key={index}
          geometry={geometry}
          onPointerMove={handlePointerMove}
          onPointerOver={handlePointerOver}
          onPointerOut={handlePointerOut}
        />
      ))}
    </>
  );
}
