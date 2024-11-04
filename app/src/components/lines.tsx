import { CatmullRomLine, Dodecahedron } from "@react-three/drei";
import { useEffect, useState } from "react";
import { socket } from "../socket";
import { useThree } from "@react-three/fiber";

import * as THREE from "three";

const findClosestPoint = (points: THREE.Vector3[], position: THREE.Vector3) => {
  const closestPoint = new THREE.Vector3();
  points.forEach((point) => {
    if (point.distanceTo(position) < closestPoint.distanceTo(position)) {
      closestPoint.copy(point);
    }
  });
  return closestPoint;
};

// TODO: ensure type consistency, every point/... should be THREE.Vector3
export const Line3D = ({
  points,
  setPoints,
  setSelectedPoint,
  isDrawing,
  colorMode,
  hoveredId, // if null, hover virtual canvas -> close line
  setIsDrawing,
  setLineLength,
}: {
  points: THREE.Vector3[];
  setPoints: any;
  setSelectedPoint: any;
  isDrawing: boolean;
  colorMode: string;
  hoveredId: number | null;
  setIsDrawing: any;
  setLineLength: (length: number) => void;
}) => {
  //   a virtual point is between every two points in the points array on the line
  const [virtualPoints, setVirtualPoints] = useState<THREE.Vector3[]>([]);
  const [lineColor, setLineColor] = useState("black");
  const [pointColor, setPointColor] = useState("black");
  const [virtualPointColor, setVirtualPointColor] = useState("darkcyan");
  const initalTriggerRef = useRef(true);

  useEffect(() => {
    // TODO: use bootstrap colors
    if (hoveredId == null && isDrawing) {
      setLineColor("#f01d23");
      setPointColor("#710000");
    } else if (colorMode === "light") {
      setLineColor("#454b66");
      setPointColor("#191308");
      setVirtualPointColor("#677db7");
    } else {
      setLineColor("#f5fdc6");
      setPointColor("#41521f");
      setVirtualPointColor("#a89f68");
    }
  }, [colorMode, hoveredId, isDrawing]);

  const handleClick = (event: any) => {
    if (!isDrawing) {
      setSelectedPoint(event.object.position.clone());
    } else {
      if (hoveredId != null) {
        const point = event.point.clone();
        setPoints([...points, point]);
      } else {
        setIsDrawing(false);
      }
    }
  };

  const handleVirtualClick = (index: number, event: any) => {
    // make the virtual point a real point and insert it at the correct position
    const newPoints = [...points];
    newPoints.splice(index + 1, 0, event.object.position.clone());
    setPoints(newPoints);
    setSelectedPoint(event.object.position.clone());
  };

  useEffect(() => {
    if (points.length < 2) {
      return;
    }
    // TODO: do not compute the curve twice
    // TODO: clean up types, reuse vector objects here
    const curve = new THREE.CatmullRomCurve3(points);

    setLineLength(curve.getLength());

    const linePoints = curve.getPoints(points.length * 20);
    const position = new THREE.Vector3();
    let _newPoints: THREE.Vector3[] = [];
    for (let i = 0; i < points.length - 1; i++) {
      position.copy(points[i]);
      position.lerp(new THREE.Vector3(...points[i + 1]), 0.5);

      _newPoints = [..._newPoints, findClosestPoint(linePoints, position)];
    }
    setVirtualPoints(_newPoints);
  }, [points]);

  useEffect(() => {
    if (initalTriggerRef.current) {
      initalTriggerRef.current = false;
      return;
    }
    if (points.length > 0) {
      // add the moving point when going from not drawing -> drawing
      // this removes a point when triggered initially
      // This should not trigger initially, so the initialTriggeRef is
      // a strange workaround
      if (isDrawing) {
        setPoints([...points, points[points.length - 1]]);
      } else {
        setPoints(points.slice(0, points.length - 1));
      }
    }
  }, [isDrawing]);

  return (
    <>
      {points.map((point, index) => (
        <Dodecahedron
          key={index}
          args={[0.1, 0]}
          position={new THREE.Vector3(...point)}
          material-color={pointColor}
          onClick={handleClick}
        />
      ))}
      {points.length >= 2 && (
        <>
          <CatmullRomLine
            points={points.map((point) => new THREE.Vector3(...point))}
            color={lineColor}
            lineWidth={2}
            segments={parseInt(points.length * 20)}
          />
          {virtualPoints.map((point, index) => (
            <Dodecahedron
              key={index}
              args={[0.1, 0]}
              position={new THREE.Vector3(...point)}
              material-color={virtualPointColor}
              onClick={(event) => handleVirtualClick(index, event)}
            />
          ))}
        </>
      )}
    </>
  );
};

import { useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { Plane } from "@react-three/drei";

export const VirtualCanvas = ({
  isDrawing,
  setPoints,
  points,
  hoveredId,
  setHoveredId,
}: {
  isDrawing: boolean;
  setPoints: any;
  points: THREE.Vector3[];
  hoveredId: number | null;
  setHoveredId: (id: number | null) => void;
}) => {
  const { camera, size } = useThree();
  const [distance, setDistance] = useState(10);
  const [canvasVisible, setCanvasVisible] = useState(false);
  const canvasRef = useRef();

  // setDistance to camera <-> last point distance

  // useEffect(() => {
  //   if (canvasRef.current) {
  //     // Set the initial size of the plane
  //     updatePlaneSize();
  //   }
  // }, [camera, size]);

  const updatePlaneSize = () => {
    const vFOV = THREE.MathUtils.degToRad(camera.fov); // Convert vertical FOV to radians
    const height = 2 * Math.tan(vFOV / 2) * distance; // Visible height
    const width = height * camera.aspect; // Visible width

    canvasRef.current.scale.set(width, height, 1); // Update the scale of the plane
  };

  const onHover = (event: any) => {
    if (isDrawing && event.object.visible) {
      if (!canvasRef.current) {
        return;
      }
      // this feature is temporarily disabled
      // if (event.shiftKey) {
      // setHoveredId(canvasRef.current);
      // // set opacity of the virtual canvas
      // setCanvasVisible(true);
      // } else {
      //   setHoveredId(null);
      //   console.log("virtual canvas");
      //   setCanvasVisible(false);
      // }

      // find the index of the closest visible point from the camera
      // if nothing is being hovered, this is the virtual canvas
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

  useFrame(() => {
    if (!isDrawing) {
      return;
    }
    if (canvasRef.current) {
      updatePlaneSize();
      // if nothing is hovered, the canvas should be visible
      canvasRef.current.visible =
        hoveredId == null || hoveredId == canvasRef.current;
    }

    if (points.length >= 2) {
      const lastPoint = points[points.length - 2];
      // the lastPoint is the one we are currently drawing
      const dist = camera.position.distanceTo(lastPoint);
      setDistance(dist);
    }

    if (canvasRef.current) {
      const direction = new THREE.Vector3(0, 0, -1).applyQuaternion(
        camera.quaternion,
      );
      const position = direction.multiplyScalar(distance).add(camera.position);
      canvasRef.current.position.copy(position);

      canvasRef.current.lookAt(camera.position);
    }
  });

  return (
    <>
      {isDrawing && (
        <Plane
          ref={canvasRef}
          args={[
            1, 1, 50, 50,
          ]} /* Initial size is [1, 1], will be scaled dynamically */
          onPointerMove={onHover}
        >
          <meshBasicMaterial
            attach="material"
            color="black"
            opacity={canvasVisible ? 0.5 : 0}
            transparent={true}
            side={THREE.DoubleSide}
          />
        </Plane>
      )}
    </>
  );
};
