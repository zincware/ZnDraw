import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import * as THREE from "three";

function createControls(camera, canvas, config, scene) {
  const controls = new OrbitControls(camera, canvas);

  // damping and auto rotation require
  // the controls to be updated each frame

  // this.controls.autoRotate = true;
  controls.enableDamping = false;

  controls.tick = function (delta) {
    if (config.selected.length > 0 && config.pressed_keys.c == true) {
      // iterate through selected and compute center

      let items = [];
      config.selected.forEach((item) => {
        items.push(scene.getObjectByName(item));
      });

      const center = items
        .reduce((a, b) => a.add(b.position), new THREE.Vector3())
        .divideScalar(items.length);
      controls.target.copy(center);
    }
    controls.update();
  };

  return controls;
}

export { createControls };
