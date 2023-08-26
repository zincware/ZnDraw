import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { TransformControls } from 'three/examples/jsm/controls/TransformControls.js';
import * as THREE from 'three';

function createControls(camera, canvas) {
  const controls = new OrbitControls(camera, canvas);

  // damping and auto rotation require
  // the controls to be updated each frame

  // this.controls.autoRotate = true;
  controls.enableDamping = false;

  controls.tick = function () {
    // if (config.selected.length > 0 && config.pressed_keys.c == true) {
    //   controls.target.copy(scene.getObjectByName("particlesGroup").get_center());
    // }
    controls.update();
  };

  return controls;
}

function createTransformControls(camera, canvas, orbit) {
  const controls = new TransformControls(camera, canvas);

  controls.addEventListener('dragging-changed', (event) => {
    orbit.enabled = !event.value;
  });
  controls.addEventListener('objectChange', () => {
    controls.object.update();
  });

  return controls;
}

export { createControls, createTransformControls };
