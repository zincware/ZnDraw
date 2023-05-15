import { Clock } from "three";

const clock = new Clock();

class Loop {
  constructor(camera, scene, renderer, stream, config) {
    this.camera = camera;
    this.scene = scene;
    this.renderer = renderer;
    this.stream = stream;
    this.updatables = [];
    this.constraint_updatables = [];
    this.config = config;
  }

  start() {
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

  tick() {
    // only call the getDelta function once per frame!
    const delta = clock.getElapsedTime();
    // console.log(
    //   `The last frame rendered in ${delta * 1000} milliseconds`,
    // );

    for (const object of this.updatables) {
      object.tick(delta);
    }

    if (delta > 1 / this.config.config.max_fps) {
      for (const object of this.constraint_updatables) {
        object.tick(this.stream.get_next_frame());
      }
      clock.start();
    }
  }
}

export { Loop };
