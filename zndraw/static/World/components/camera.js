import * as THREE from "three";

function createCamera() {
  // const camera = new THREE.PerspectiveCamera(
  //   35, // fov = Field Of View
  //   1, // aspect ratio (dummy value)
  //   0.1, // near clipping plane
  //   1000, // far clipping plane
  // );
  let width = window.innerWidth;
  let height = window.innerHeight;

  const camera = new THREE.OrthographicCamera(
    width / -2,
    width / 2,
    height / 2,
    height / -2,
  );

  // move the camera back so we can view the scene
  camera.position.set(0, 0, 20);
  const cameraLight = new THREE.PointLight(0xffffff, 1.0);
  camera.add(cameraLight);

  return camera;
}

export { createCamera };
