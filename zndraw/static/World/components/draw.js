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
    this.name = "canvas3D";

    this.add(plane, wireframe);
  };
}

export class Line3D extends THREE.Group {
  constructor() {
    super();
    this.anchorPoints = [];

    this.ARC_SEGMENTS = 200;

    const geometry = new THREE.BufferGeometry();
    const material = new THREE.LineBasicMaterial({ color: 0xff0000 });
    this.line = new THREE.Line(geometry, material);
    this.curve = undefined;

    this.control_points = new THREE.Group();
    this.add(this.line, this.control_points);

    // call updateControlPoints() when x is pressed
    document.addEventListener("keydown", (event) => {
      if (event.key === "x") {
        this.updateControlPoints();
      }
    });
  }

  addPoint(position) {
    // const geometry = new THREE.SphereGeometry(1.4, 32, 32);
    // const material = new THREE.MeshBasicMaterial({ color: "#000000" });
    // const sphere = new THREE.Mesh(geometry, material);
    // sphere.position.copy(position);

    this.anchorPoints.push(position);
    this.updateLine();
    // this.add(sphere);
  }

  updateLine() {
    if (this.anchorPoints.length < 2) {
      return;
    }
    this.curve = new THREE.CatmullRomCurve3(this.anchorPoints);

    const points = this.curve.getPoints(this.ARC_SEGMENTS);
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    this.line.geometry = geometry;
    this.updateControlPoints();
  }

  movePointer(position) {
    this.anchorPoints[this.anchorPoints.length - 1].copy(position);
    this.updateLine();
  }

  updateControlPoints() {
    // remove all children from three group control_points
    // while (this.control_points.children.length > 0) {
    //   this.control_points.remove(this.control_points.children[0]);
    // }
    // create 5 spheres along the full distance of the curve

    // create a control point between each pair of anchor points
    // the point has to be on the curve

    // for (let i = 0; i < this.anchorPoints.length - 1; i++) {
    //   const anchorPoint1 = this.anchorPoints[i];
    //   const anchorPoint2 = this.anchorPoints[i + 1];
    //   const controlPoint = new THREE.Vector3().addVectors(anchorPoint1, anchorPoint2).divideScalar(2);

    //   const geometry = new THREE.SphereGeometry(0.2, 32, 32);
    //   const material = new THREE.MeshBasicMaterial({ color: "#000000" });
    //   const sphere = new THREE.Mesh(geometry, material);
    //   sphere.position.copy(controlPoint);
    //   this.control_points.add(sphere);
    // }

    // const n_points = Math.floor(this.curve.getLength());

    // this.curve.getSpacedPoints(n_points).forEach((point) => {
    //   const geometry = new THREE.SphereGeometry(0.2, 32, 32);
    //   const material = new THREE.MeshBasicMaterial({ color: "#000000" });
    //   const sphere = new THREE.Mesh(geometry, material);
    //   sphere.position.copy(point);
    //   this.control_points.add(sphere);
    // });
  }
}
