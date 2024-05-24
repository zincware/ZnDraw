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
    const old_step = this.step;
    this.step = step;
    let succeeded = [];

    for (const object of this.step_updatables) {
      const success = object.step(step);
      if (success) {
        // console.log("success", object);
        succeeded.push(object);
      } else {
        // console.log("failed", object);
        break;
      }
    }
    if (succeeded.length == this.step_updatables.length) {
      return true;
    } else {
      // undo changes if not all succeeded
      // invesitgate if this is needed
      for (const object of succeeded) {
        object.step(step - 1);
      }
      this.step = old_step;
      return false;
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
