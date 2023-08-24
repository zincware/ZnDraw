import * as THREE from "three";
import { TransformControls } from "three/examples/jsm/controls/TransformControls.js";

class Selection {
  constructor(camera, scene, socket, line3D, renderer, controls) {
    this.camera = camera;
    this.scene = scene;
    this.socket = socket;
    this.controls = controls;

    this.line3D = line3D;

    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();
    this.selection = [];
    this.transform_controls = new TransformControls(
      camera,
      renderer.domElement,
    );

    // I don't like this here! Add in world.
    this.scene.add(this.transform_controls);

    const onPointerMove = this.onPointerMove.bind(this);

    // event on backspace
    document.addEventListener("keydown", (event) => {
      if (event.key === "Backspace") {
        this.line3D.removePointer(this.transform_controls.object);
        this.transform_controls.detach();
      }
      if (event.key === "Escape") {
        this.transform_controls.detach();
        if (this._drawing) {
          this._drawing = false;
          this.line3D.removePointer();
          window.removeEventListener("pointermove", onPointerMove);
        }
      }
    });

    this.transform_controls.addEventListener(
      "dragging-changed",
      function (event) {
        controls.enabled = !event.value;
      },
    );
    this.transform_controls.addEventListener("objectChange", () => {
      // mesh -> anchorPoints -> Line3D
      this.transform_controls.object.parent.parent.updateLine();
    });

    this._drawing = false;

    window.addEventListener("pointerdown", this.onPointerDown.bind(this));

    // use x keypress to toggle the attachment of onPointerMove
    document.addEventListener("keydown", (event) => {
      if (event.key === "x") {
        if (this._drawing) {
          this._drawing = false;
          this.line3D.removePointer();
          window.removeEventListener("pointermove", onPointerMove);
        } else {
          this._drawing = true;
          this.line3D.addPointer();
          this.transform_controls.detach();
          window.addEventListener("pointermove", onPointerMove);
        }
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

    return this.raycaster.intersectObjects(this.scene.children, true);
  }

  /**
   * Drawing raycaster
   */
  onPointerMove(event) {
    const intersects = this.getIntersections();
    const particlesGroup = this.scene.getObjectByName("particlesGroup");
    for (let i = 0; i < intersects.length; i++) {
      const object = intersects[i].object;
      if (object.name === "canvas3D") {
        this.line3D.movePointer(intersects[i].point.clone());
        break;
      }
      if (particlesGroup.children.includes(object.parent)) {
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
      if (this._drawing) {
        if (object.name === "canvas3D") {
          this.line3D.addPoint(intersects[i].point.clone());
          break;
        }
        if (particlesGroup.children.includes(object.parent)) {
          this.line3D.addPoint(intersects[i].point.clone());
          break;
        }
      } else {
        if (object.parent.name === "AnchorPoints") {
          this.transform_controls.attach(object);
        } else {
          if (particlesGroup.children.includes(object.parent)) {
            if (this.selection.includes(object.parent.name)) {
              this.selection = this.selection.filter(
                (x) => x !== object.parent.name,
              );
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
  }
}

export { Selection };
