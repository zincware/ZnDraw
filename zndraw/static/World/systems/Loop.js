import { Clock } from 'three';

const clock = new Clock();
const constraint_clock = new Clock();

class Loop {
  constructor(camera, scene, renderer, renderer2d, socket) {
    this.camera = camera;
    this.scene = scene;
    this.renderer = renderer;
    this.renderer2d = renderer2d;
    this.tick_updatables = [];
    this.step_updatables = [];

    this.step = 0;

    this.config = { max_fps: 60 };
    // Update the config object when the server sends a new one
    socket.on('config', (data) => {
      this.config = data;
      console.log('Loop config:');
      console.log(this.config);
    });
  }

  start() {
    this.renderer.setAnimationLoop(() => {
      this.tick();
    });
    this.setStep(this.step);
  }

  stop() {
    this.renderer.setAnimationLoop(null);
  }

  setStep(step) {
    this.step = step;

    for (const object of this.step_updatables) {
      object.step(step);
    }
  }

  tick() {
    this.renderer.render(this.scene, this.camera);
    this.renderer2d.render(this.scene, this.camera);

    for (const object of this.tick_updatables) {
      object.tick();
    }
  }
}

export { Loop };
