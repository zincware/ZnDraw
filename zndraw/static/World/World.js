import { createCamera } from './components/camera.js';
import { createLights } from './components/lights.js';
import { createScene } from './components/scene.js';

import { createControls, createTransformControls } from './systems/controls.js';
import { createRenderer, create2DRenderer } from './systems/renderer.js';
import { Resizer } from './systems/Resizer.js';
import { Loop } from './systems/Loop.js';
import { ParticlesGroup } from './components/particles.js';
import { Selection } from './systems/select.js';

import { Line3D, Canvas3D } from './components/draw.js';

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

class Player {
  constructor(world, cache) {
    this.world = world;
    this.playing = false;
    this.fps = 30;
    this.cache = cache;

    // toggle playing on spacebar
    document.addEventListener('keydown', (event) => {
      if (document.activeElement === document.body && event.code === 'Space') {
        this.toggle();
      }
      if (
        document.activeElement === document.body
        && event.code === 'ArrowRight'
      ) {
        this.go_forward();
      }
      // on arrow left go backward
      if (
        document.activeElement === document.body
        && event.code === 'ArrowLeft'
      ) {
        this.go_backward();
      }
      // on arrow up go forward 10 % of the length
      if (
        document.activeElement === document.body
        && event.code === 'ArrowUp'
      ) {
        this.go_forward(parseInt(this.cache.get_length() / 10));
      }
      // on arrow down go backward 10 % of the length
      if (
        document.activeElement === document.body
        && event.code === 'ArrowDown'
      ) {
        this.go_backward(parseInt(this.cache.get_length() / 10));
      }
    });

    const slider = document.getElementById('frame-slider');
    slider.focus();

    slider.oninput = function () {
      document.getElementById(
        'info',
      ).innerHTML = `${slider.value} / ${slider.max}`;
      world.setStep(this.value);
    };
  }

  toggle() {
    this.playing = !this.playing;
    if (this.playing) this.play();
  }

  go_forward(step = 1) {
    let new_step = this.world.getStep() + step;
    if (new_step >= this.cache.get_length()) {
      new_step = step - 1;
    }
    this.world.setStep(new_step);
  }

  go_backward(step = 1) {
    let new_step = this.world.getStep() - step;
    if (new_step < 0) {
      new_step = this.cache.get_length() - step;
    }
    this.world.setStep(new_step);
  }

  play() {
    if (this.playing) {
      this.go_forward();
      setTimeout(() => this.play(), 1000 / this.fps);
    }
  }
}

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
    this.line3D = new Line3D(camera, renderer);
    const canvas3D = new Canvas3D();

    selection = new Selection(
      camera,
      scene,
      socket,
      this.line3D,
      renderer,
      controls,
    );

    const light = createLights();

    scene.add(particles, light, camera, this.line3D, canvas3D); // index, transform_controls

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
    step = parseInt(step);
    loop.setStep(step);
    const slider = document.getElementById('frame-slider');
    slider.value = step;
    document.getElementById(
      'info',
    ).innerHTML = `${slider.value} / ${slider.max}`;
  }

  getStep() {
    return loop.step;
  }

  getSelection() {
    return selection.selection;
  }

  getPoints() {
    // create a list of positions from this.line3D.anchorPoints;
    return this.line3D.anchorPoints.children.map((x) => x.position);
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

export { World, Player };
