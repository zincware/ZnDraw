import * as THREE from "three";
import { TransformControls } from "three/examples/jsm/controls/TransformControls.js";
import { centerCamera, duplicateAnchorPoints } from "./events.js";

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

    // add a function to the line3D callbacks that is triggered if the line is changed
    // why is this in selection?
    this.line3D.onLineChange = (emit = true) => {
      // let points = this.line3D.anchorPoints.children.map((x) => x.position);
      // // convert x, y, z to [x, y, z]
      // points = points.map((x) => [x.x, x.y, x.z]);
      // this.socket.emit("room:points:set", { 0: points });
      if (emit) {
        //
        // TODO: this can't work because it fixes the points in place
        // should only be done on click or not saved to redis.
      }
      // if transform controls is not attached to a point, detach it
      if (
        this.transform_controls.object &&
        this.transform_controls.object.name === "AnchorPoint"
      ) {
        if (
          !this.line3D.anchorPoints.children.includes(
            this.transform_controls.object,
          )
        ) {
          this.transform_controls.detach();
        }
      }
    };

    this.raycaster = new THREE.Raycaster();
    this.pointer = new THREE.Vector2();
    this.selection = [];
    this.transform_controls = new TransformControls(
      camera,
      renderer.domElement,
    );
    this.transform_controls.size = 0.75;

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
        this.line3D.onLineChange();
      }
    });

    window.addEventListener("click", this.onClick.bind(this));
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
      const params = document.getElementById(
        "selection-json-editor",
      ).parameters;
      this.socket.emit("selection:run", params);
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
      if (!this.line3D.pointer) {
        this.line3D.pointer = this.line3D.addPointer();
      }
      const position = particleIntersects[0].point.clone();
      this.line3D.movePointer(position, event.clientX, event.clientY);
      this.line3D.changeLineColor(0x000000);
      this.line3D.changeLastPointColor(0x000000);
    } else if (
      canvasIntersects.length > 0 &&
      canvasIntersects[0].object.name === "canvas3D"
    ) {
      if (!this.line3D.pointer) {
        this.line3D.pointer = this.line3D.addPointer();
      }
      const position = canvasIntersects[0].point.clone();
      this.line3D.movePointer(position, event.clientX, event.clientY);
      this.line3D.changeLineColor(0x000000);
      this.line3D.changeLastPointColor(0x000000);
    } else {
      if (!this.line3D.pointer) {
        this.line3D.pointer = this.line3D.addPointer();
      }
      // if (this.line3D.pointer) {
      // add a plane with the position of the pointer perpendicular to the camera
      const plane = new THREE.Plane();
      plane.setFromNormalAndCoplanarPoint(
        this.camera.getWorldDirection(plane.normal),
        this.line3D.pointer.position,
      );
      const position = new THREE.Vector3();
      this.raycaster.ray.intersectPlane(plane, position);
      this.line3D.movePointer(position, event.clientX, event.clientY);
      this.line3D.changeLineColor(0xff0000);
      this.line3D.changeLastPointColor(0xff0000);
      // }
    }
    return false;
  }

  onClick(event) {
    const elements = document.elementsFromPoint(event.clientX, event.clientY);
    // if neither first or second element is scene-container, then it's a UI element
    // first one can be the canvas
    if (
      elements[0].id !== "scene-container" &&
      elements[1].id !== "scene-container"
    ) {
      return;
    }

    // detect double click
    if (event.detail === 2) {
      this.onDoubleClick(event);
      return;
    }
    const particlesGroup = this.scene.getObjectByName("particlesGroup");
    const anchorPoints = this.scene.getObjectByName("AnchorPoints");
    const canvas3D = this.scene.getObjectByName("canvas3D");
    const virtualPoints = this.scene.getObjectByName("VirtualPoints");
    // const canvasIntersects = this.getIntersections(this
    const anchorPointsIntersects = this.getIntersections(anchorPoints);
    const particleIntersects = this.getIntersections(particlesGroup);
    const canvasIntersects = this.getIntersections(canvas3D);
    const virtualPointsIntersects = this.getIntersections(virtualPoints);

    if (virtualPointsIntersects.length > 0) {
      const position = virtualPointsIntersects[0].point.clone();
      const point = this.line3D.addPoint(
        position,
        virtualPointsIntersects[0].object.index + 1,
      );
      this.transform_controls.attach(point);
      if (this._drawing) {
        document.getElementById("drawingSwitch").click();
      }
    } else if (this._drawing) {
      if (particleIntersects.length > 0) {
        const position = particleIntersects[0].point.clone();
        this.line3D.pointer = this.line3D.addPoint(position);
      } else if (
        canvasIntersects.length > 0 &&
        canvasIntersects[0].object.name === "canvas3D"
      ) {
        const position = canvasIntersects[0].point.clone();
        this.line3D.pointer = this.line3D.addPoint(position);
      } else {
        // exit drawing mode
        document.getElementById("drawingSwitch").click();
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
        this.socket.emit("room:selection:set", { 0: particlesGroup.selection });
      }
    }
    const selectionUpdate = new CustomEvent("selection:update", {
      detail: { selection: particlesGroup.selection },
    });
    document.dispatchEvent(selectionUpdate);

    // update the points
    let points = this.line3D.anchorPoints.children.map((x) => x.position);
    // convert x, y, z to [x, y, z]
    points = points.map((x) => [x.x, x.y, x.z]);
    this.socket.emit("room:points:set", { 0: points });
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

    // check if document.getElementById("drawingSwitch") is checked and update this._drawing
    document.getElementById("drawingSwitch").addEventListener("change", () => {
      this._drawing = document.getElementById("drawingSwitch").checked;
      if (this._drawing) {
        window.addEventListener("pointermove", onPointerMove);
      } else {
        if (this.line3D.pointer) {
          this.line3D.removePointer(this.line3D.pointer);
          this.line3D.pointer = undefined;
        }
        window.removeEventListener("pointermove", onPointerMove);
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
            if (this.line3D.pointer) {
              this.line3D.removePointer(this.line3D.pointer);
              this.line3D.pointer = undefined;
            }
            window.removeEventListener("pointermove", onPointerMove);
            document.getElementById("drawingSwitch").checked = false;
          } else {
            this._drawing = true;
            this.transform_controls.detach();
            window.addEventListener("pointermove", onPointerMove);
            document.getElementById("drawingSwitch").checked = true;
          }
        }

        // make a copy of the currently selected point when d is pressed
        if (event.key === "d") {
          duplicateAnchorPoints.bind(this)();
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
            this.socket.emit("modifier:run", {
              method: { discriminator: "Delete" },
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
            if (this.line3D.pointer) {
              this.line3D.removePointer(this.line3D.pointer);
              this.line3D.pointer = undefined;
            }
            window.removeEventListener("pointermove", onPointerMove);
            document.getElementById("drawingSwitch").checked = false;
          }
        }
      }
    });

    const onPointerMove = this.onPointerMove.bind(this);
  }
}

export { Selection };
