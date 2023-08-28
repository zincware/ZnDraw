import { createCamera } from './components/camera.js';
import { createLights } from './components/lights.js';
import { createScene } from './components/scene.js';

import { createControls, createTransformControls } from './systems/controls.js';
import { createRenderer, create2DRenderer } from './systems/renderer.js';
import { Resizer } from './systems/Resizer.js';
import { Loop } from './systems/Loop.js';
import { ParticlesGroup, ParticleIndexGroup } from './components/particles.js';
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
  constructor(world, cache, socket) {
    this.world = world;
    this.playing = false;
    this.fps = 30;
    this.cache = cache;
    this.loop = false;

    socket.on('view:set', (index) => {
      this.world.setStep(index);
    });

    socket.on('view:play', () => {
      this.playing = true;
    });

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

  setLoop(loop) {
    this.loop = loop;
  }

  toggle() {
    if ((!this.playing) && this.world.getStep() == this.cache.get_length() - 1) {
      this.world.setStep(0);
    }
    this.playing = !this.playing;
  }

  go_forward(step = 1) {
    let new_step = this.world.getStep() + step;
    if (new_step >= this.cache.get_length()) {
      if (this.loop) {
        new_step = step - 1;
      } else {
        new_step = this.world.getStep();
        this.playing = false;
      }
    }
    this.world.setStep(new_step);
  }

  go_backward(step = 1) {
    let new_step = this.world.getStep() - step;
    if (new_step < 0) {
      if (this.loop) {
        new_step = this.cache.get_length() - step;
      } else {
        new_step = 0;
      }
    }
    this.world.setStep(new_step);
  }

  tick() {
    if (this.playing) {
      this.go_forward();
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

    this.player = new Player(this, cache, socket);

    container.append(renderer.domElement);

    // const transform_controls = createTransformControls(
    //   camera,
    //   renderer.domElement,
    //   controls,
    // );

    this.particles = new ParticlesGroup(socket, cache);
    this.line3D = new Line3D(camera, renderer);
    const canvas3D = new Canvas3D();
    const index_grp = new ParticleIndexGroup(this.particles);

    this.selection = new Selection(
      camera,
      scene,
      socket,
      this.line3D,
      renderer,
      controls,
      cache,
      this,
    );

    const light = createLights();

    scene.add(this.particles, light, camera, this.line3D, canvas3D, index_grp); // index, transform_controls

    // attach the canvas3D to the camera while t is pressed. attach to the scene when released
    document.addEventListener('keydown', (event) => {
      if (event.key == 't') {
        if (camera.children.includes(canvas3D)) {
          scene.attach(canvas3D);
          document.getElementById('alertBoxDrawing').style.display = 'none';
        } else {
          camera.attach(canvas3D);
          document.getElementById('alertBoxDrawing').style.display = 'block';
        }
      }
    });

    loop.tick_updatables.push(controls, this.player);
    loop.step_updatables.push(this.particles, this.selection, index_grp);

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
  rebuild(resolution, material, wireframe, simulation_box, bonds) {
    this.particles.rebuild(resolution, material, wireframe, simulation_box, bonds);
    this.setStep(loop.step);
  }

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
    return this.selection.selection;
  }

  getLineData() {
    // create a list of positions from this.line3D.anchorPoints;
    const points = this.line3D.anchorPoints.children.map((x) => x.position);
    let segments = [];
    try {
      segments = this.line3D.curve.getSpacedPoints(this.line3D.ARC_SEGMENTS).map((x) => x.toArray());
    } catch (error) {
      // console.log(error);
    }

    return { points, segments };
  }
}

export { World };
