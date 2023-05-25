import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { TransformControls } from "three/examples/jsm/controls/TransformControls.js";
import * as THREE from "three";

function createControls(camera, canvas, config, scene) {
  const controls = new OrbitControls(camera, canvas);

  // damping and auto rotation require
  // the controls to be updated each frame

  // this.controls.autoRotate = true;
  controls.enableDamping = false;

  controls.tick = function (delta) {
    if (config.selected.length > 0 && config.pressed_keys.c == true) {
      controls.target.copy(scene.getObjectByName("particleGroup").get_center());
    }
    controls.update();
  };

  return controls;
}

function createTransformControls(camera, canvas, orbit) {
  const controls = new TransformControls(camera, canvas);

  controls.addEventListener("dragging-changed", function (event) {
    orbit.enabled = !event.value;
  });
  controls.addEventListener("objectChange", function () {
    controls.object.update();
  });

  return controls;
}

export { createControls, createTransformControls };
