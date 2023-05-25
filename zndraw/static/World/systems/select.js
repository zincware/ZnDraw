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

    // iterate itersections until we find a particle
    for (let i = 0; i < intersects.length; i++) {
      const particleGroup = this.scene.getObjectByName("particleGroup");
      const object = intersects[i].object;

      if (object.name == "drawCanvas") {
        object.click(intersects[i].point);
      }

      if (particleGroup.children.includes(object.parent)) {
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
