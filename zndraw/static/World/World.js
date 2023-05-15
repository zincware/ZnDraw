import { createCamera } from "./components/camera.js";
import { createLights } from "./components/lights.js";
import { createScene } from "./components/scene.js";

import { createControls } from "./systems/controls.js";
import { createRenderer } from "./systems/renderer.js";
import { Resizer } from "./systems/Resizer.js";
import { Loop } from "./systems/Loop.js";
import { Stream } from "./systems/Stream.js";
import { createParticleGroup } from "./components/particles.js";

// These variables are module-scoped: we cannot access them
// from outside the module
let camera;
let renderer;
let scene;
let loop;
let stream;

class World {
  constructor(container) {
    camera = createCamera();
    scene = createScene();
    renderer = createRenderer();
    stream = new Stream();
    loop = new Loop(camera, scene, renderer, stream);

    container.append(renderer.domElement);

    const controls = createControls(camera, renderer.domElement);

    const particles = createParticleGroup();
    const light = createLights();

    scene.add(particles, light);

    // disable mesh rotation
    loop.constraint_updatables.push(particles);
    loop.updatables.push(controls);

    const resizer = new Resizer(container, camera, renderer);
  }

  render() {
    // draw a single frame
    renderer.render(scene, camera);
  }
  start() {
    loop.start();
  }

  stop() {
    loop.stop();
  }
}

export { World };
