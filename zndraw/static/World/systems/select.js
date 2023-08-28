import * as THREE from 'three';
import { TransformControls } from 'three/examples/jsm/controls/TransformControls.js';

let scroll_timer = null;

class Selection {
  constructor(camera, scene, socket, line3D, renderer, controls, cache, world) {
    this.camera = camera;
    this.scene = scene;
    this.socket = socket;
    this.controls = controls;
    this.cache = cache;
    this.world = world;

    this.shift_pressed = false;

    document.addEventListener('keydown', (event) => {
      if (event.key === 'Shift') {
        this.shift_pressed = true;
      }
    });
    document.addEventListener('keyup', (event) => {
      if (event.key === 'Shift') {
        this.shift_pressed = false;
      }
    });

    this.controls.getCenter = () => {
      const particlesGroup = this.scene.getObjectByName('particlesGroup');
      return particlesGroup.get_center();
    };

    this.socket.on('selection:run', (data) => {
      this.selection = data;
      this.step();
    });

    window.addEventListener('wheel', this.onWheel.bind(this));

    // use c keypress to center the camera on the selection
    document.addEventListener('keydown', (event) => {
      if (event.key === 'c') {
        if (this.controls.enablePan) {
          document.getElementById('alertBoxCamera').style.display = 'block';

          // get the first object that is selected
          const particlesGroup = this.scene.getObjectByName('particlesGroup');

          particlesGroup.children.every((x) => {
            if (this.selection.includes(x.name)) {
              this.controls.target = x.position;
              this.controls.enablePan = false;
              return false;
              // TODO: don't use the first but the COM of the selection
            }
            return true;
          });
        } else {
          document.getElementById('alertBoxCamera').style.display = 'none';
          this.controls.target = this.controls.target.clone();
          this.controls.enablePan = true;
        }
      }
    });
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
    document.addEventListener('keydown', (event) => {
      if (event.key === 'Backspace' && document.activeElement === document.body) {
        this.line3D.removePointer(this.transform_controls.object);
        this.transform_controls.detach();
      }
      if (event.key === 'Escape' && document.activeElement === document.body) {
        this.transform_controls.detach();
        if (this._drawing) {
          this._drawing = false;
          this.line3D.removePointer();
          window.removeEventListener('pointermove', onPointerMove);
        }
      }
    });

    this.transform_controls.addEventListener('dragging-changed', (event) => {
      controls.enabled = !event.value;
    });
    this.transform_controls.addEventListener('objectChange', () => {
      if (this.transform_controls.object.name === 'Canvas3DGroup') {
        // nothing to do
      } else {
        // mesh -> anchorPoints -> Line3D
        this.transform_controls.object.parent.parent.updateLine();
      }
    });

    this._drawing = false;

    window.addEventListener('pointerdown', this.onPointerDown.bind(this));
    window.addEventListener('dblclick', this.onDoubleClick.bind(this));

    // use x keypress to toggle the attachment of onPointerMove
    document.addEventListener('keydown', (event) => {
      if (event.key === 'x') {
        if (this._drawing) {
          this._drawing = false;
          this.line3D.removePointer();
          window.removeEventListener('pointermove', onPointerMove);
        } else {
          this._drawing = true;
          this.line3D.addPointer();
          this.transform_controls.detach();
          window.addEventListener('pointermove', onPointerMove);
        }
      }
    });
  }

  step() {
    const particlesGroup = this.scene.getObjectByName('particlesGroup');
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

  onDoubleClick(event) {
    const intersects = this.getIntersections();
    const particlesGroup = this.scene.getObjectByName('particlesGroup');
    const selection = [];
    for (let i = 0; i < intersects.length; i++) {
      const { object } = intersects[i];
      if (particlesGroup.children.includes(object.parent)) {
        selection.push(object.parent.name);

        const selectionOptions = document.getElementById('selection-select');

        this.socket.emit('selection:run', {
          name: selectionOptions.options[selectionOptions.selectedIndex].text,
          params: {},
          atoms: this.cache.get(this.world.getStep()),
          selection,
        });
        break;
      }
    }
  }

  /**
   * Drawing raycaster
   */
  onPointerMove(event) {
    const intersects = this.getIntersections();
    const particlesGroup = this.scene.getObjectByName('particlesGroup');
    for (let i = 0; i < intersects.length; i++) {
      const { object } = intersects[i];
      if (object.name === 'canvas3D') {
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
      const particlesGroup = this.scene.getObjectByName('particlesGroup');
      const { object } = intersects[i];
      if (this._drawing) {
        if (object.name === 'canvas3D') {
          this.line3D.addPoint(intersects[i].point.clone());
          break;
        }
        if (particlesGroup.children.includes(object.parent)) {
          this.line3D.addPoint(intersects[i].point.clone());
          break;
        }
      } else if (object.parent.name === 'AnchorPoints') {
        this.transform_controls.attach(object);
      } else if (particlesGroup.children.includes(object.parent)) {
        if (this.shift_pressed) {
          if (this.selection.includes(object.parent.name)) {
            this.selection = this.selection.filter(
              (x) => x !== object.parent.name,
            );
            object.parent.set_selection(false);
          } else {
            this.selection.push(object.parent.name);
            object.parent.set_selection(true);
          }
        } else {
          if (this.selection.includes(object.parent.name)) {
            this.selection = [];
          } else {
            this.selection = [object.parent.name];
          }
          object.parent.set_selection(true);
          this.step();
        }
        break; // only (de)select one particle
      }
    }
  }

  onWheel(event) {
    const intersections = this.getIntersections();
    for (let i = 0; i < intersections.length; i++) {
      const { object } = intersections[i];
      if (object.name === 'canvas3D') {
        console.log('scrolled on canvas3D');
        // there must be a better way to disable scrolling while over the canvas
        this.controls.enabled = false;
        if (scroll_timer) {
          clearTimeout(scroll_timer);
        }
        scroll_timer = setTimeout(() => {
          this.controls.enabled = true;
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

export { Selection };
