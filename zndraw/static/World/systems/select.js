import * as THREE from "three";

class Selection {
  constructor(camera, scene, socket) {
    this.camera = camera;
    this.scene = scene;
    this.socket = socket;

    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();
    this.selection = [];

    window.addEventListener("pointerdown", this.onPointerDown.bind(this));
  }

  step() {
    const particlesGroup = this.scene.getObjectByName("particlesGroup");
    // iterate through all children ids that are in the selection and update them
    particlesGroup.children.forEach((x) => {
      if (this.selection.includes(x.name)) {
        x.set_selection(true);
      } else {
        x.set_selection(false);
      }
    });
  }

  onPointerDown(event) {
    this.pointer.x = (event.clientX / window.innerWidth) * 2 - 1;
    this.pointer.y = -(event.clientY / window.innerHeight) * 2 + 1;

    this.raycaster.setFromCamera(this.pointer, this.camera);

    const intersects = this.raycaster.intersectObjects(
      this.scene.children,
      true,
    );

    // iterate intersections until we find a particle
    for (let i = 0; i < intersects.length; i++) {
      const particlesGroup = this.scene.getObjectByName("particlesGroup");
      const object = intersects[i].object;
      if (particlesGroup.children.includes(object.parent)) {
        if (this.selection.includes(object.parent.name)) {
          this.selection = this.selection.filter((x) => x !== object.parent.name);
          object.parent.set_selection(false);
        } else {
          this.selection.push(object.parent.name);
          object.parent.set_selection(true);
        }
        this.socket.emit("selection", this.selection);
        break; // only (de)select one particle
      }
    }
  }
}

export { Selection };
