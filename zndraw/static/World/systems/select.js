import * as THREE from "three";

class Selection {
  constructor(camera, scene, config) {
    this.camera = camera;
    this.scene = scene;
    this.config = config;

    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();

    window.addEventListener("pointerdown", this.onPointerDown.bind(this));
  }

  async onPointerDown(event) {
    this.pointer.x = (event.clientX / window.innerWidth) * 2 - 1;
    this.pointer.y = -(event.clientY / window.innerHeight) * 2 + 1;

    this.raycaster.setFromCamera(this.pointer, this.camera);

    const intersects = this.raycaster.intersectObjects(
      this.scene.children,
      true,
    );

    if (intersects.length > 0) {
      const object = intersects[0].object;
      if (object.parent.parent.name != "particleGroup") {
        return;
      }
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
    }
  }
}

export { Selection };
