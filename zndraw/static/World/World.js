import { createCamera } from "./components/camera.js";
import { createLights } from "./components/lights.js";
import { createScene } from "./components/scene.js";

import { createControls } from "./systems/controls.js";
import { createRenderer, create2DRenderer } from "./systems/renderer.js";
import { Resizer } from "./systems/Resizer.js";
import { Loop } from "./systems/Loop.js";
import {
  ParticlesGroup,
  ParticleIndexGroup,
  CellGroup,
} from "./components/particles.js";
import { Selection } from "./systems/select.js";

import { Line3D, Canvas3D } from "./components/draw.js";
import { centerCamera } from "./systems/events.js";

// These variables are module-scoped: we cannot access them
// from outside the module
let camera;
let renderer;
let renderer2d;
let scene;
let loop;
let controls;

class Player {
  constructor(world, cache, socket) {
    this.world = world;
    this.playing = false;
    this.fps = 60;
    this.cache = cache;
    this.loop = false;

    socket.on("scene:set", (index) => {
      this.world.setStep(index);
    });

    socket.on("view:play", () => {
      this.playing = true;
      this.play();
    });

    // detect playBtn click
    document.getElementById("playBtn").addEventListener("click", () => {
      this.toggle();
      // if playing set the icon to <i class="fa-solid fa-pause"></i> else to <i class="fa-solid fa-play"></i>
      if (this.playing) {
        document.getElementById("playBtn").innerHTML =
          '<i class="fa-solid fa-pause"></i>';
      }
      if (!this.playing) {
        document.getElementById("playBtn").innerHTML =
          '<i class="fa-solid fa-play"></i>';
      }
    });
    // detect forwardBtn click
    document.getElementById("forwardBtn").addEventListener("click", () => {
      this.go_forward();
    });
    // detect backwardBtn click
    document.getElementById("backwardBtn").addEventListener("click", () => {
      this.go_backward();
    });

    // toggle playing on spacebar
    document.addEventListener("keydown", (event) => {
      if (document.activeElement === document.body && event.code === "Space") {
        this.toggle();
      }
      if (
        document.activeElement === document.body &&
        event.code === "ArrowRight"
      ) {
        this.go_forward();
      }
      // on arrow left go backward
      if (
        document.activeElement === document.body &&
        event.code === "ArrowLeft"
      ) {
        this.go_backward();
      }
      // on arrow up go forward 10 % of the length
      if (
        document.activeElement === document.body &&
        event.code === "ArrowUp"
      ) {
        this.go_forward(parseInt(this.cache.get_length() / 10));
      }
      // on arrow down go backward 10 % of the length
      if (
        document.activeElement === document.body &&
        event.code === "ArrowDown"
      ) {
        this.go_backward(parseInt(this.cache.get_length() / 10));
      }
    });

    const slider = document.getElementById("frame-slider");
    slider.focus();

    slider.oninput = function () {
      document.getElementById(
        "info",
      ).innerHTML = `${slider.value} / ${slider.max}`;
      world.setStep(this.value);
    };
  }

  setLoop(loop) {
    this.loop = loop;
  }

  toggle() {
    if (!this.playing && this.world.getStep() == this.cache.get_length() - 1) {
      this.world.setStep(0);
    }
    this.playing = !this.playing;
    if (this.playing) this.play();
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
    cache.attachWorld(this);

    loop = new Loop(camera, scene, renderer, renderer2d);
    controls = createControls(camera, renderer.domElement);

    this.player = new Player(this, cache, socket);

    container.append(renderer.domElement);

    this.particles = new ParticlesGroup(socket, cache);
    this.line3D = new Line3D(camera, renderer);
    this.index_grp = new ParticleIndexGroup(this.particles, camera);

    this.cell_grp = new CellGroup(cache);

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

    const canvas3D = new Canvas3D(this.particles);

    const light = createLights();

    scene.add(
      this.particles,
      light,
      camera,
      this.line3D,
      canvas3D,
      this.index_grp,
      this.cell_grp,
    ); // index, transform_controls

    // attach the canvas3D to the camera while t is pressed. attach to the scene when released
    document.addEventListener("keydown", (event) => {
      if (document.activeElement === document.body) {
        if (event.key == "f") {
          if (camera.children.includes(canvas3D)) {
            scene.attach(canvas3D);
            document.getElementById("alertBoxDrawing").style.display = "none";
          } else {
            camera.attach(canvas3D);
            document.getElementById("alertBoxDrawing").style.display = "block";
          }
        }
      }
    });

    loop.tick_updatables.push(controls, this.index_grp);
    loop.step_updatables.push(this.particles, this.index_grp, this.cell_grp);

    const resizer = new Resizer(container, camera, renderer, renderer2d);

    this.step = loop.step;
    this.socket = socket;

    this.socket.on(
      "scene:points",
      function (callback) {
        const { points, segments } = this.getLineData();
        callback(points);
      }.bind(this),
    );

    this.socket.on(
      "scene:segments",
      function (callback) {
        const { points, segments } = this.getLineData();
        callback(segments);
      }.bind(this),
    );

    this.socket.on(
      "scene:step",
      function (callback) {
        callback(this.getStep());
      }.bind(this),
    );

    this.socket.on(
      "selection:get",
      function (callback) {
        callback(this.getSelection());
      }.bind(this),
    );
    // renderer.render(scene, camera);
  }

  /**
   * Start the event loop
   */
  start() {
    loop.start();
    // center camera if there are particles
    const particles = this.particles.cache.get(0);
    if (!(particles === undefined || particles.length === 0)) {
      centerCamera(this.selection.controls, this.particles);
    }
  }

  /**
   * Rebuild all objects in the scene
   */
  rebuild(
    resolution,
    material,
    wireframe,
    simulation_box,
    bonds,
    label_offset,
    particle_size,
    bonds_size,
    fps,
  ) {
    this.particles.rebuild(
      resolution,
      material,
      wireframe,
      bonds,
      particle_size,
      bonds_size,
    );
    this.cell_grp.set_visibility(simulation_box);
    this.setStep(loop.step);
    this.index_grp.rebuild(label_offset);
    this.player.fps = fps;
  }

  setStep(step) {
    step = parseInt(step);
    loop.setStep(step);
    const slider = document.getElementById("frame-slider");
    slider.value = step;
    document.getElementById(
      "info",
    ).innerHTML = `${slider.value} / ${slider.max}`;
  }

  getStep() {
    return loop.step;
  }

  getSelection() {
    return this.particles.selection;
  }

  getLineData() {
    // create a list of positions from this.line3D.anchorPoints;
    const points = this.line3D.anchorPoints.children.map((x) => x.position);
    let segments = [];
    try {
      if (points.length > 1) {
        segments = this.line3D.curve
          .getSpacedPoints(this.line3D.ARC_SEGMENTS)
          .map((x) => x.toArray());
      }
    } catch (error) {
      // console.log(error);
    }

    return { points, segments };
  }
}

export { World };
