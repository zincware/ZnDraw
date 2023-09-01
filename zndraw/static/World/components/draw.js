import * as THREE from "three";

export class Canvas3D extends THREE.Group {
  constructor(particlesGroup) {
    super();
    this.name = "Canvas3DGroup";

    let material;

    let geometry;

    const drawAddCanvasBtn = document.getElementById("drawAddCanvas");

    drawAddCanvasBtn.addEventListener("click", () => {
      this.removeCanvas();

      material = new THREE.MeshBasicMaterial({
        color: drawAddCanvasBtn.parameters.color,
        side: THREE.DoubleSide,
        transparent: true,
        opacity: drawAddCanvasBtn.parameters.opacity,
      });

      // TODO use json-forms to create dynamic forms for each geometry

      // const geometry_name = document.getElementById('drawCanvasSelect').value;
      const params = drawAddCanvasBtn.parameters.geometry;
      const geometry_name = params.method;

      if (geometry_name === "PlaneGeometry") {
        geometry = new THREE.PlaneGeometry(params.width, params.height);
      } else if (geometry_name === "BoxGeometry") {
        geometry = new THREE.BoxGeometry(
          params.width,
          params.height,
          params.depth,
        );
      } else if (geometry_name === "CircleGeometry") {
        geometry = new THREE.CircleGeometry(params.radius, 32);
      } else if (geometry_name === "ConeGeometry") {
        geometry = new THREE.ConeGeometry(params.radius, params.height, 32);
      } else if (geometry_name === "CylinderGeometry") {
        geometry = new THREE.CylinderGeometry(
          params.radius_top,
          params.radius_bottom,
          params.height,
          32,
        );
      } else if (geometry_name === "DodecahedronGeometry") {
        geometry = new THREE.DodecahedronGeometry(params.radius);
      } else if (geometry_name === "IcosahedronGeometry") {
        geometry = new THREE.IcosahedronGeometry(params.radius);
      } else if (geometry_name === "OctahedronGeometry") {
        geometry = new THREE.OctahedronGeometry(params.radius);
      } else if (geometry_name === "RingGeometry") {
        geometry = new THREE.RingGeometry(
          params.inner_radius,
          params.outer_radius,
          32,
        );
      } else if (geometry_name === "SphereGeometry") {
        geometry = new THREE.SphereGeometry(params.radius, 32, 32);
      } else if (geometry_name === "TetrahedronGeometry") {
        geometry = new THREE.TetrahedronGeometry(params.radius);
      } else if (geometry_name === "TorusGeometry") {
        geometry = new THREE.TorusGeometry(params.radius, params.tube, 32, 100);
      } else if (geometry_name === "TorusKnotGeometry") {
        geometry = new THREE.TorusKnotGeometry(
          params.radius,
          params.tube,
          100,
          16,
        );
      }

      const plane = new THREE.Mesh(geometry, material);
      plane.name = "canvas3D";

      if (drawAddCanvasBtn.parameters.wireframe) {
        let wireframeGeometry;

        if (true) {
          const thresholdAngle = 1;
          wireframeGeometry = new THREE.EdgesGeometry(geometry, thresholdAngle);
        } else {
          wireframeGeometry = new THREE.WireframeGeometry(geometry);
        }
        const wireframeMaterial = new THREE.LineBasicMaterial({
          color: 0x000000,
          transparent: true,
          opacity: 0.5,
        });
        const wireframe = new THREE.LineSegments(
          wireframeGeometry,
          wireframeMaterial,
        );
        wireframe.name = "canvas3D-wireframe";
        this.add(wireframe);
      }
      this.add(plane);
      if (particlesGroup.selection.length > 0) {
        const matrix = new THREE.Matrix4();
        const dummy = new THREE.Object3D();
        particlesGroup.particles_mesh.getMatrixAt(
          particlesGroup.selection[0],
          matrix,
        );
        matrix.decompose(dummy.position, dummy.quaternion, dummy.scale);
        this.position.copy(dummy.position);
      }
    });

    document
      .getElementById("drawRemoveCanvas")
      .addEventListener("click", () => {
        this.removeCanvas();
      });
  }

  removeCanvas() {
    this.remove(this.getObjectByName("canvas3D"));
    this.remove(this.getObjectByName("canvas3D-wireframe"));
  }
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

  addPoint(position, index) {
    const geometry = new THREE.IcosahedronGeometry(0.1, 0);
    const material = new THREE.MeshPhongMaterial({
      color: "#000000",
      shininess: 100,
    });
    const sphere = new THREE.Mesh(geometry, material);
    sphere.name = "AnchorPoint";
    sphere.position.copy(position);
    if (index !== undefined) {
      const allMeshes = [...this.anchorPoints.children];
      this.anchorPoints.clear();

      allMeshes.splice(index, 0, sphere);
      this.anchorPoints.add(...allMeshes);
    } else {
      this.anchorPoints.add(sphere);
    }

    this.updateLine();

    this.ARC_SEGMENTS = this.anchorPoints.children.length * 20;

    return sphere;
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
