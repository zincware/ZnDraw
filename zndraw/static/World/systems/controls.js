import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";

function createControls(camera, canvas) {
  const controls = new OrbitControls(camera, canvas);

  controls.enableDamping = false;

  controls.tick = function () {
    controls.update();
  };

  return controls;
}

export { createControls };
