class Loop {
  constructor(camera, scene, renderer, renderer2d) {
    this.camera = camera;
    this.scene = scene;
    this.renderer = renderer;
    this.renderer2d = renderer2d;
    this.tick_updatables = [];
    this.step_updatables = [];

    this.step = 0;
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
