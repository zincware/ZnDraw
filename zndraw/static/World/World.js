import { createCamera } from "./components/camera.js";
import { createLights } from "./components/lights.js";
import { createScene } from "./components/scene.js";

import { createControls, createTransformControls } from "./systems/controls.js";
import { createRenderer, create2DRenderer } from "./systems/renderer.js";
import { Resizer } from "./systems/Resizer.js";
import { Loop } from "./systems/Loop.js";
import { Stream } from "./systems/Stream.js";
import {
  ParticleGroup,
  createParticleGroup,
  createIndexGroup,
} from "./components/particles.js";
import { Selection } from "./systems/select.js";

import { Curve3D } from "./components/draw.js";

// These variables are module-scoped: we cannot access them
// from outside the module
let camera;
let renderer;
let renderer2d;
let scene;
let loop;
let controls
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

    const particles = new ParticleGroup(socket, cache);

    const light = createLights();

    scene.add(particles, light, camera); // index, transform_controls

    loop.updatables.push(controls);
    loop.constraint_updatables.push(particles);

    const resizer = new Resizer(container, camera, renderer, renderer2d);

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
  rebuild() {

  }

  setStep(step) {
    loop.setStep(step);
  }

}
    

// class World {
//   constructor(container, cache) {
   
//     selection = new Selection(camera, scene, config);

//     
//     const transform_controls = createTransformControls(
//       camera,
//       renderer.domElement,
//       controls,
//     );

//     
//     const index = createIndexGroup(particles);
//     
//     const curve = new Curve3D(scene, transform_controls, config, particles);

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
