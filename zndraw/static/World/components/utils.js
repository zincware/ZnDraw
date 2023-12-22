import * as THREE from "three";

export const findClosestPoint = (points, position) => {
  const closestPoint = new THREE.Vector3();
  points.forEach((point) => {
    if (point.distanceTo(position) < closestPoint.distanceTo(position)) {
      closestPoint.copy(point);
    }
  });
  return closestPoint;
};
