import { Clock } from "three";

const clock = new Clock();
const constraint_clock = new Clock();

class Loop {
  constructor(camera, scene, renderer, renderer2d, socket) {
    this.camera = camera;
    this.scene = scene;
    this.renderer = renderer;
    this.renderer2d = renderer2d;
    this.updatables = [];
    this.constraint_updatables = [];

    this.step = 0;

    this.config = { "max_fps": 60 };
    // Update the config object when the server sends a new one
    socket.on("config", (data) => {
      this.config = data;
      console.log("Loop config:");
      console.log(this.config);
    });
  }

  start() {
    this.tick();
    this.setStep(this.step);

    this.renderer.setAnimationLoop(() => {
      // tell every animated object to tick forward one frame
      this.tick();

      // render a frame
      this.renderer.render(this.scene, this.camera);
    });
  }

  stop() {
    this.renderer.setAnimationLoop(null);
  }

  setStep(step) {
    this.step = step;
    for (const object of this.constraint_updatables) {
      object.step(step);
    }
  }

  tick() {
    // only call the getDelta function once per frame!
    // split into a tick() and a frame() function Maybe trigger frame() via socket?
    this.renderer2d.render(this.scene, this.camera);

    // TODO only update necessary objects. Don't split into updatables and constraint_updatables
    for (const object of this.updatables) {
      object.tick();
    }

    // if (clock.getElapsedTime() > 1 / this.config.max_fps) {
    //   for (const object of this.constraint_updatables) {
    //     object.step(0);  // This somehow must update the current frame to display
    //   }
    //   clock.start();
    // }
  }
}

export { Loop };
