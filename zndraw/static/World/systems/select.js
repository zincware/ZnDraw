import * as THREE from "three";

class Selection {
  constructor(camera, scene, config, controls) {
    this.camera = camera;
    this.scene = scene;
    this.config = config;
    this.controls = controls;

    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();

    this._disable_controls_timeout_id = null;
    this._wheel_target = null;

    window.addEventListener("pointerdown", this.onPointerDown.bind(this));
    window.addEventListener("wheel", this.onWheel.bind(this));
  }

  async getIntersected(event) {
    this.pointer.x = (event.clientX / window.innerWidth) * 2 - 1;
    this.pointer.y = -(event.clientY / window.innerHeight) * 2 + 1;

    this.raycaster.setFromCamera(this.pointer, this.camera);

    const intersects = this.raycaster.intersectObjects(
      this.scene.children,
      true,
    );

    return intersects;
  }

  async onWheel(event) {
    const intersects = await this.getIntersected(event);

    function renewTimeOut() {
      if (typeof this._disable_controls_timeout_id === "number") {
        clearTimeout(this._disable_controls_timeout_id);
      }
      this._disable_controls_timeout_id = setTimeout(() => {
        this.controls.enabled = true;
        this._wheel_target = null;
      }, 500);
    }

    if (this._wheel_target === null) {
      for (let i = 0; i < intersects.length; i++) {
        const object = intersects[i].object;

        if (object.name == "drawCanvas") {
          this.controls.enabled = false;
          this._wheel_target = object;

          renewTimeOut.bind(this)();

          console.log(event.deltaY);
          object.scale.set(
            object.scale.x + event.deltaY * 0.0005,
            object.scale.y + event.deltaY * 0.0005,
            object.scale.z + event.deltaY * 0.0005,
          );
          break;
        }
      }
    } else {
      renewTimeOut.bind(this)();
      this._wheel_target.scale.set(
        this._wheel_target.scale.x + event.deltaY * 0.0005,
        this._wheel_target.scale.y + event.deltaY * 0.0005,
        this._wheel_target.scale.z + event.deltaY * 0.0005,
      );
    }
  }

  async onPointerDown(event) {
    const intersects = await this.getIntersected(event);

    // iterate itersections until we find a particle
    for (let i = 0; i < intersects.length; i++) {
      const particleGroup = this.scene.getObjectByName("particleGroup");
      const object = intersects[i].object;

      if (object.name == "drawCanvas") {
        object.click(intersects[i].point);
        break;
      }

      if (
        object.userData.type === "particle" ||
        object.userData.type === "bond"
      ) {
        object.parent.click();
        await fetch("select", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            selected_ids: this.config.selected,
            step: this.config.step,
            method: document.getElementById("selection-method").value,
          }),
        })
          .then((response) => response.json())
          .then((data) => (this.config.selected = data["selected_ids"]));
        break;
      }
    }
  }
}

export { Selection };
