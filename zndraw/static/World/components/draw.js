import * as THREE from 'three';

export class Canvas3D extends THREE.Group {
  constructor() {
    super();
    this.name = 'Canvas3DGroup';
    const material = new THREE.MeshBasicMaterial({
      color: '#cccccc',
      side: THREE.DoubleSide,
      transparent: true,
      opacity: 0.5,
    });

    let geometry;

    document.getElementById('drawAddCanvas').addEventListener('click', () => {
      this.remove(this.getObjectByName('canvas3D'));
      this.remove(this.getObjectByName('canvas3D-wireframe'));

      // TODO use json-forms to create dynamic forms for each geometry

      const geometry_name = document.getElementById('drawCanvasSelect').value;
      if (geometry_name === 'PlaneGeometry') {
        geometry = new THREE.PlaneGeometry(10, 10);
      } else if (geometry_name === 'BoxGeometry') {
        geometry = new THREE.BoxGeometry(10, 10, 10);
      } else if (geometry_name === 'CircleGeometry') {
        geometry = new THREE.CircleGeometry(5, 32);
      } else if (geometry_name === 'ConeGeometry') {
        geometry = new THREE.ConeGeometry(5, 20, 32);
      } else if (geometry_name === 'CylinderGeometry') {
        geometry = new THREE.CylinderGeometry(5, 5, 20, 32);
      } else if (geometry_name === 'DodecahedronGeometry') {
        geometry = new THREE.DodecahedronGeometry(5);
      } else if (geometry_name === 'IcosahedronGeometry') {
        geometry = new THREE.IcosahedronGeometry(5);
      } else if (geometry_name === 'OctahedronGeometry') {
        geometry = new THREE.OctahedronGeometry(5);
      } else if (geometry_name === 'RingGeometry') {
        geometry = new THREE.RingGeometry(1, 5, 32);
      } else if (geometry_name === 'SphereGeometry') {
        geometry = new THREE.SphereGeometry(5, 32, 32);
      } else if (geometry_name === 'TetrahedronGeometry') {
        geometry = new THREE.TetrahedronGeometry(5);
      } else if (geometry_name === 'TorusGeometry') {
        geometry = new THREE.TorusGeometry(5, 1, 32, 100);
      } else if (geometry_name === 'TorusKnotGeometry') {
        geometry = new THREE.TorusKnotGeometry(5, 1, 100, 16);
      }

      const wireframeGeometry = new THREE.WireframeGeometry(geometry);
      const wireframeMaterial = new THREE.LineBasicMaterial({ color: 0x000000 });
      const wireframe = new THREE.LineSegments(
        wireframeGeometry,
        wireframeMaterial,
      );
      const plane = new THREE.Mesh(geometry, material);
      plane.name = 'canvas3D';
      wireframe.name = 'canvas3D-wireframe';

      this.add(plane, wireframe);
    });

    document.getElementById('drawRemoveCanvas').addEventListener('click', () => {
      // clean previous canvas
      this.remove(this.getObjectByName('canvas3D'));
      this.remove(this.getObjectByName('canvas3D-wireframe'));
    });
  }
}

export class Line3D extends THREE.Group {
  constructor() {
    super();
    this.anchorPoints = new THREE.Group();
    this.anchorPoints.name = 'AnchorPoints';

    this.ARC_SEGMENTS = 200;

    const geometry = new THREE.BufferGeometry();
    const material = new THREE.LineBasicMaterial({ color: 0x000000 });
    this.line = new THREE.Line(geometry, material);
    this.curve = undefined;

    this.add(this.line, this.anchorPoints);
  }

  addPoint(position) {
    const geometry = new THREE.IcosahedronGeometry(0.1, 0);
    const material = new THREE.MeshPhongMaterial({
      color: '#000000',
      shininess: 100,
    });
    const sphere = new THREE.Mesh(geometry, material);
    sphere.position.copy(position);
    this.anchorPoints.add(sphere);

    this.updateLine();

    this.ARC_SEGMENTS = this.anchorPoints.children.length * 20;
  }

  removePointer(object) {
    // remove last anchor point
    if (object) {
      this.anchorPoints.remove(object);
    } else {
      this.anchorPoints.remove(
        this.anchorPoints.children[this.anchorPoints.children.length - 1],
      );
    }

    this.updateLine();
  }

  addPointer() {
    this.addPoint(new THREE.Vector3(0, 0, 0));
  }

  updateLine() {
    if (this.anchorPoints.children.length < 2) {
      // remove the line
      this.line.geometry = new THREE.BufferGeometry();
      return;
    }

    // get list of positions from anchor points
    this.curve = new THREE.CatmullRomCurve3(
      this.anchorPoints.children.map((x) => x.position),
    );

    const points = this.curve.getPoints(this.ARC_SEGMENTS);
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    this.line.geometry = geometry;
  }

  movePointer(position) {
    // this.anchorPoints[this.anchorPoints.length - 1].copy(position);
    this.anchorPoints.children[
      this.anchorPoints.children.length - 1
    ].position.copy(position);
    // create red material
    this.updateLine();
  }
}
