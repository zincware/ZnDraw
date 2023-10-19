import * as THREE from "three";
import { TransformControls } from "three/examples/jsm/controls/TransformControls.js";
import { centerCamera } from "./events.js";

let scroll_timer = null;

class Selection {
  constructor(camera, scene, socket, line3D, renderer, controls, cache, world) {
    this.camera = camera;
    this.scene = scene;
    this.socket = socket;
    this.controls = controls;
    this.cache = cache;
    this.world = world;
    this._drawing = false;

    this.line3D = line3D;

    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();
    this.selection = [];
    this.transform_controls = new TransformControls(
      camera,
      renderer.domElement,
    );

    // change mode of transform_controls when pressing t, iterate between translate, rotate, scale
    document.addEventListener("keydown", (event) => {
      if (document.activeElement === document.body) {
        if (event.key === "t") {
          if (this.transform_controls.mode === "translate") {
            this.transform_controls.setMode("rotate");
          } else if (this.transform_controls.mode === "rotate") {
            this.transform_controls.setMode("scale");
          } else if (this.transform_controls.mode === "scale") {
            this.transform_controls.setMode("translate");
          }
        }
        if (event.key === "f") {
          this.transform_controls.detach();
        }
      }
    });

    // I don't like this here! Add in world.
    this.scene.add(this.transform_controls);

    this.shift_pressed = false;
    this.ctrl_pressed = false;

    this._setupKeyboardEvents();

    this.controls.getCenter = (selection) => {
      const particlesGroup = this.scene.getObjectByName("particlesGroup");
      return particlesGroup.get_center(selection);
    };

    this.socket.on("selection:set", (data) => {
      const particlesGroup = this.scene.getObjectByName("particlesGroup");
      particlesGroup.selection = data;
      particlesGroup.step();
    });

    window.addEventListener("wheel", this.onWheel.bind(this));

    this.transform_controls.addEventListener("dragging-changed", (event) => {
      controls.enabled = !event.value;
    });
    this.transform_controls.addEventListener("objectChange", () => {
      if (this.transform_controls.object.name === "Canvas3DGroup") {
        // nothing to do
      } else {
        // mesh -> anchorPoints -> Line3D
        this.transform_controls.object.parent.parent.updateLine();
      }
    });

    window.addEventListener("pointerdown", this.onPointerDown.bind(this));
    window.addEventListener("dblclick", this.onDoubleClick.bind(this));
  }

  getIntersections(object) {
    this.pointer.x = (event.clientX / window.innerWidth) * 2 - 1;
    this.pointer.y = -(event.clientY / window.innerHeight) * 2 + 1;

    this.raycaster.setFromCamera(this.pointer, this.camera);
    if (object) {
      return this.raycaster.intersectObject(object, true);
    } else {
      return this.raycaster.intersectObjects(this.scene.children, true);
    }
  }

  onDoubleClick(event) {
    const particlesGroup = this.scene.getObjectByName("particlesGroup");

    const particleIntersects = this.getIntersections(particlesGroup);
    if (particleIntersects.length > 0) {
      const instanceId = particleIntersects[0].instanceId;
      particlesGroup.click(
        instanceId,
        this.shift_pressed,
        particleIntersects[0].object,
      );
      const params = document.getElementById(
        "selection-json-editor",
      ).parameters;
      console.log(new Date().toISOString(), "running selection");
      this.socket.emit("selection:run", {
        params: params,
      });
    }
  }

  /**
   * Drawing raycaster
   */
  onPointerMove(event) {
    const particlesGroup = this.scene.getObjectByName("particlesGroup");
    const canvas3D = this.scene.getObjectByName("canvas3D");

    const particleIntersects = this.getIntersections(particlesGroup);
    const canvasIntersects = this.getIntersections(canvas3D);

    if (particleIntersects.length > 0) {
      const position = particleIntersects[0].point.clone();
      this.line3D.movePointer(position);
    } else if (canvasIntersects.length > 0) {
      if (canvasIntersects[0].object.name === "canvas3D") {
        const position = canvasIntersects[0].point.clone();
        this.line3D.movePointer(position);
      }
    }
    // } else {
    //   this.line3D.removePointer();
    // }
    return false;
  }

  onPointerDown(event) {
    const particlesGroup = this.scene.getObjectByName("particlesGroup");
    const anchorPoints = this.scene.getObjectByName("AnchorPoints");
    const canvas3D = this.scene.getObjectByName("canvas3D");
    // const canvasIntersects = this.getIntersections(this
    const anchorPointsIntersects = this.getIntersections(anchorPoints);
    const particleIntersects = this.getIntersections(particlesGroup);
    const canvasIntersects = this.getIntersections(canvas3D);

    if (this._drawing) {
      if (particleIntersects.length > 0) {
        const position = particleIntersects[0].point.clone();
        this.line3D.addPoint(position);
      } else if (canvasIntersects.length > 0) {
        if (canvasIntersects[0].object.name === "canvas3D") {
          const position = canvasIntersects[0].point.clone();
          this.line3D.addPoint(position);
        }
      }
    } else {
      if (anchorPointsIntersects.length > 0) {
        const object = anchorPointsIntersects[0].object;
        if (object.name === "AnchorPoint") {
          this.transform_controls.attach(object);
        }
      } else if (particleIntersects.length > 0) {
        const instanceId = particleIntersects[0].instanceId;
        particlesGroup.click(
          instanceId,
          this.shift_pressed,
          particleIntersects[0].object,
        );
      }
    }
  }

  onWheel(event) {
    if (this.shift_pressed) {
      const intersections = this.getIntersections();
      for (let i = 0; i < intersections.length; i++) {
        const { object } = intersections[i];
        if (object.name === "canvas3D") {
          // there must be a better way to disable scrolling while over the canvas
          this.controls.enableZoom = false;
          if (scroll_timer) {
            clearTimeout(scroll_timer);
          }
          scroll_timer = setTimeout(() => {
            this.controls.enableZoom = true;
          }, 500);

          this.transform_controls.attach(object.parent);

          object.parent.children.forEach((x) => {
            x.scale.set(
              x.scale.x + event.deltaY * 0.0005,
              x.scale.y + event.deltaY * 0.0005,
              x.scale.z + event.deltaY * 0.0005,
            );
          });

          break;
        }
      }
    }
  }

  _setupKeyboardEvents() {
    document.addEventListener("keydown", (event) => {
      if (event.key === "Shift") {
        this.shift_pressed = true;
      }
      if (event.key === "Control") {
        this.ctrl_pressed = true;
      }
    });
    document.addEventListener("keyup", (event) => {
      if (event.key === "Shift") {
        this.shift_pressed = false;
      }
      if (event.key === "Control") {
        this.ctrl_pressed = false;
      }
    });

    // use c keypress to center the camera on the selection
    document.addEventListener("keydown", (event) => {
      if (document.activeElement === document.body) {
        const particlesGroup = this.scene.getObjectByName("particlesGroup");
        if (event.key === "c") {
          centerCamera(this.controls, particlesGroup);
        }
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

        // make a copy of the currently selected point when d is pressed
        if (event.key === "d") {
          // TODO shift added point a bit + insert at correct position
          const transform_object = this.transform_controls.object;
          if (transform_object.name === "AnchorPoint" && !this._drawing) {
            const index =
              this.line3D.anchorPoints.children.indexOf(transform_object);

            let new_pos;
            if (index > 0) {
              // Add the point between the current and the previous point
              const obj_before = this.line3D.anchorPoints.children[index - 1];
              new_pos = obj_before.position
                .clone()
                .sub(transform_object.position)
                .multiplyScalar(0.5)
                .add(transform_object.position);
            } else {
              // No previous point, add the point at the same position
              new_pos = transform_object.position.clone();
            }

            const point = this.line3D.addPoint(new_pos, index);
            this.transform_controls.detach();
            this.transform_controls.attach(point);
          }
        }

        if (event.key === "Backspace") {
          // remove pointer if transform_controls is attached to it
          if (this.transform_controls.object) {
            if (this.transform_controls.object.name === "AnchorPoint") {
              this.line3D.removePointer(this.transform_controls.object);
              this.transform_controls.detach();
            } else if (
              this.transform_controls.object.name === "Canvas3DGroup"
            ) {
              this.transform_controls.object.removeCanvas();
              this.transform_controls.detach();
            }
          } else if (particlesGroup.selection.length > 0) {
            const { points, segments } = this.world.getLineData();
            console.log(new Date().toISOString(), "running modifier");
            this.socket.emit("modifier:run", {
              params: { method: { method: "Delete" } },
              url: window.location.href,
            });
            // particlesGroup.click();
          } else {
            this.line3D.removePointer();
          }
        }
        if (event.key === "Escape") {
          this.transform_controls.detach();
          if (this._drawing) {
            this._drawing = false;
            this.line3D.removePointer();
            window.removeEventListener("pointermove", onPointerMove);
          }
        }
      }
    });

    const onPointerMove = this.onPointerMove.bind(this);
  }
}

export { Selection };
