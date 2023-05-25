import * as THREE from "three";

export class Canvas3D {
  constructor(scene, transformControls, config, particleGroup, camera, curve) {
    this.scene = scene;
    this.transformControls = transformControls;
    this.config = config;
    this.particleGroup = particleGroup;
    this.camera = camera;
    this.curve = curve;

    document.getElementById("drawAddCanvas").onclick = () => {
      this.createCanvas();
    }
    document.getElementById("drawRemoveCanvas").onclick = () => {
      this.removeCanvas();
    }

    this.canvas = null;
  }

  click(point) {
    console.log(point);
    this.curve.createAnchorPoint(point);
  }

  createCanvas() {
    const geometry = new THREE.PlaneGeometry(10, 10);
    const material = new THREE.MeshBasicMaterial({ color: 0xcccccc, side: THREE.DoubleSide, transparent: true, opacity: 0.5 });
    const wireframeGeometry = new THREE.WireframeGeometry(geometry);
    const wireframeMaterial = new THREE.LineBasicMaterial({ color: 0x000000});
    const wireframe = new THREE.LineSegments(wireframeGeometry, wireframeMaterial);
    const plane = new THREE.Mesh(geometry, material);
    plane.add(wireframe);
    plane.name = "drawCanvas";
    plane.position.copy(this.particleGroup.get_center());
    plane.lookAt(this.camera.position);
    plane.click = this.click.bind(this);
    this.canvas = plane;
    this.scene.add(plane);
  }

  removeCanvas() {
    if (this.canvas) {
      this.scene.remove(this.canvas);
      this.canvas = null;
    }
  }

}

export class Curve3D {
  constructor(scene, transformControls, config, particleGroup) {
    this.scene = scene;
    this.transformControls = transformControls;
    this.config = config;
    this.particleGroup = particleGroup;

    this.ARC_SEGMENTS = 200;
    this.curve;

    this.createCurve();
    this.anchorPoints = new THREE.Group();

    this.scene.add(this.curve);
    this.scene.add(this.anchorPoints);

    document.getElementById("drawAddAnchor").onclick = () => {
      this.createAnchorPoint();
    };
    document.getElementById("drawRemoveLine").onclick = () => {
      this.removeCurve();
    };

    document.getElementById("drawDetach").onclick = () => {
      this.transformControls.detach();
    };
  }

  createCurve() {
    // const curve = new THREE.CatmullRomCurve3( this.points );

    // const points = curve.getPoints( this.ARC_SEGMENTS );
    const geometry = new THREE.BufferGeometry();

    const material = new THREE.LineBasicMaterial({ color: 0xff0000 });

    // Create the final object to add to the scene
    const curveObject = new THREE.Line(geometry, material);

    this.curve = curveObject;
  }

  createAnchorPoint(position = undefined) {
    if (position === undefined) {
      const position = this.particleGroup.get_center();
    }

    const geometry = new THREE.SphereGeometry(0.2, 32, 32);
    const material = new THREE.MeshBasicMaterial({ color: "#000000" });
    const sphere = new THREE.Mesh(geometry, material);
    sphere.position.copy(position);

    this.anchorPoints.add(sphere);
    this.transformControls.attach(sphere);

    // save the points here for passing to flask
    this.config.draw_vectors = Array.from(
      this.anchorPoints.children,
      (child) => child.position,
    );

    sphere.update = () => {
      if (this.config.draw_vectors.length < 2) {
        return;
      }
      this.curve.geometry.dispose();
      const curve = new THREE.CatmullRomCurve3(this.config.draw_vectors);

      const points = curve.getPoints(this.ARC_SEGMENTS);
      const geometry = new THREE.BufferGeometry().setFromPoints(points);
      this.curve.geometry = geometry;
    };

    sphere.update();
  }

  removeCurve() {
    while (this.anchorPoints.children.length > 0) {
      this.anchorPoints.remove(this.anchorPoints.children[0]);
    }
    this.transformControls.detach();
    this.curve.geometry.dispose();
    this.curve.geometry = new THREE.BufferGeometry();
    this.config.draw_vectors = [];
  }
}
