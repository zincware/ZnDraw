import { createCamera } from "./components/camera.js";
import { createLights } from "./components/lights.js";
import { createScene } from "./components/scene.js";

import { createControls, createTransformControls } from "./systems/controls.js";
import { createRenderer, create2DRenderer } from "./systems/renderer.js";
import { Resizer } from "./systems/Resizer.js";
import { Loop } from "./systems/Loop.js";
import { Stream } from "./systems/Stream.js";
import { ParticlesGroup } from "./components/particles.js";
import { Selection } from "./systems/select.js";

import { Line3D, Canvas3D } from "./components/draw.js";

// These variables are module-scoped: we cannot access them
// from outside the module
let camera;
let renderer;
let renderer2d;
let scene;
let loop;
let controls;
let cache;
let selection;

class World {
  constructor(container, cache, socket) {
    camera = createCamera();
    scene = createScene();
    renderer = createRenderer();
    renderer2d = create2DRenderer();

    cache = cache;

    loop = new Loop(camera, scene, renderer, renderer2d, socket);
    controls = createControls(camera, renderer.domElement);

    container.append(renderer.domElement);

    // const transform_controls = createTransformControls(
    //   camera,
    //   renderer.domElement,
    //   controls,
    // );

    const particles = new ParticlesGroup(socket, cache);
    const line3D = new Line3D(camera, renderer);
    const canvas3D = new Canvas3D();

    selection = new Selection(
      camera,
      scene,
      socket,
      line3D,
      renderer,
      controls,
    );

    const light = createLights();

    scene.add(particles, light, camera, line3D, canvas3D); // index, transform_controls

    loop.tick_updatables.push(controls);
    loop.step_updatables.push(particles, selection);

    const resizer = new Resizer(container, camera, renderer, renderer2d);

    this.step = loop.step;
    this.socket = socket;

    // renderer.render(scene, camera);
  }

  /**
   * Start the event loop
   */
  start() {
    loop.start();
  }

  /**
   * Rebuild all objects in the scene
   */
  rebuild() {}

  setStep(step) {
    loop.setStep(step);
    this.socket.emit("step", {"step": step});
  }
}

// class World {
//   constructor(container, cache) {

//     selection = new Selection(camera, scene, config);

//

//
//     const index = createIndexGroup(particles);
//
//

//     window.addEventListener("keydown", (event) => {
//       if (event.isComposing || event.key === "i") {
//         index.show();
//       }
//     });
//     // remove index group when i is released
//     window.addEventListener("keyup", (event) => {
//       if (event.isComposing || event.key === "i") {
//         index.hide();
//       }
//     });

//     document.getElementById("reset").onclick = () => {
//       this.deleteCache();
//       this.rebuild();
//     };

//

//     // disable mesh rotation
//     loop.constraint_updatables.push(index);
//

//
//   }

//   render() {
//     // draw a single frame
//
//   }
//   start() {
//
//   }

//   stop() {
//     loop.stop();
//   }

//   rebuild() {
//     // remove all children from particles
//     scene.children[0].clear();
//   }

// }

export { World };
