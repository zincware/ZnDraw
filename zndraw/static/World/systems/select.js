import * as THREE from "three";

class Selection {
  constructor(camera, scene, socket, line3D) {
    this.camera = camera;
    this.scene = scene;
    this.socket = socket;

    this.line3D = line3D;

    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();
    this.selection = [];

    window.addEventListener("pointerdown", this.onPointerDown.bind(this));

    const onPointerMove = this.onPointerMove.bind(this);
    // use x keypress to toggle the attachment of onPointerMove
    document.addEventListener("keydown", (event) => {
      if (event.key === "x") {
        window.addEventListener("pointermove", onPointerMove);
      }
      if (event.key === "c") {
        console.log("Remove pointermove event listener");
        window.removeEventListener("pointermove", onPointerMove);
      }
    });
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

  getIntersections() {
    this.pointer.x = (event.clientX / window.innerWidth) * 2 - 1;
    this.pointer.y = -(event.clientY / window.innerHeight) * 2 + 1;

    this.raycaster.setFromCamera(this.pointer, this.camera);

    return this.raycaster.intersectObjects(
      this.scene.children,
      true,
    );
  }

  /**
   * Drawing raycaster
   */
  onPointerMove(event) {
    console.log("onPointerMove");
    const intersects = this.getIntersections();
    const particlesGroup = this.scene.getObjectByName("particlesGroup");
    for (let i = 0; i < intersects.length; i++) {
      const object = intersects[i].object;
      if (object.parent.name === "canvas3D") {
        this.line3D.movePointer(intersects[i].point.clone());
        break;
      }
      if (particlesGroup.children.includes(object.parent)){
        this.line3D.movePointer(intersects[i].point.clone());
        break;
      }
    }
    return false;
  }

  onPointerDown(event) {
    const intersects = this.getIntersections();

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
          this.line3D.addPoint(object.parent.position.clone());
        }
        this.socket.emit("selection", this.selection);
        break; // only (de)select one particle
      }
    }
  }
}

export { Selection };
