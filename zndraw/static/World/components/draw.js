import * as THREE from "three";

export class Canvas3D extends THREE.Group {
  constructor() {
    super();
    const geometry = new THREE.PlaneGeometry(10, 10);
    const material = new THREE.MeshBasicMaterial({
      color: "#cccccc",
      side: THREE.DoubleSide,
      transparent: true,
      opacity: 0.5,
    });
    const wireframeGeometry = new THREE.WireframeGeometry(geometry);
    const wireframeMaterial = new THREE.LineBasicMaterial({ color: 0x000000 });
    const wireframe = new THREE.LineSegments(
      wireframeGeometry,
      wireframeMaterial,
    );
    const plane = new THREE.Mesh(geometry, material);
    plane.name = "canvas3D";

    this.add(plane, wireframe);
  };
}

export class Line3D extends THREE.Group {
  constructor() {
    super();
    this.anchorPoints = new THREE.Group();
    this.anchorPoints.name = "AnchorPoints";

    this.ARC_SEGMENTS = 200;

    const geometry = new THREE.BufferGeometry();
    const material = new THREE.LineBasicMaterial({ color: 0x000000 });
    this.line = new THREE.Line(geometry, material);
    this.curve = undefined;

    this.add(this.line, this.anchorPoints);
  }

  addPoint(position) {
    const geometry = new THREE.IcosahedronGeometry( 0.1, 0 );
    const material = new THREE.MeshPhongMaterial({ color: "#000000", shininess: 100});
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
      this.anchorPoints.remove(this.anchorPoints.children[this.anchorPoints.children.length - 1]);
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
    this.curve = new THREE.CatmullRomCurve3(this.anchorPoints.children.map((x) => x.position)); 

    const points = this.curve.getPoints(this.ARC_SEGMENTS);
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    this.line.geometry = geometry;

  }

  movePointer(position) {
    // this.anchorPoints[this.anchorPoints.length - 1].copy(position);
    this.anchorPoints.children[this.anchorPoints.children.length - 1].position.copy(position);
    // create red material
    this.updateLine();
  }
}
