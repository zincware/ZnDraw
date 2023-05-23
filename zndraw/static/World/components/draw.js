import * as THREE from 'three';


export class Curve3D {
    constructor(scene, transformControls) {
        this.scene = scene;
        this.transformControls = transformControls;

        this.ARC_SEGMENTS = 200;
        this.curve;

        this.createCurve();
        this.anchorPoints = new THREE.Group();

        this.scene.add(this.curve);
        this.scene.add(this.anchorPoints);        
    }

    createCurve() {
        // const curve = new THREE.CatmullRomCurve3( this.points );
        
        // const points = curve.getPoints( this.ARC_SEGMENTS );
        const geometry = new THREE.BufferGeometry();
        
        const material = new THREE.LineBasicMaterial( { color: 0xff0000 } );
        
        // Create the final object to add to the scene
        const curveObject = new THREE.Line( geometry, material );

        this.curve = curveObject;
    }

    createAnchorPoint() {
        const position = new THREE.Vector3(0, 0, 0);
        const geometry = new THREE.SphereGeometry(0.2, 32, 32);
        const material = new THREE.MeshBasicMaterial({ color: "#000000" });
        const sphere = new THREE.Mesh(geometry, material);
        sphere.position.copy(position);

        sphere.update = () => {
            this.curve.geometry.dispose();
            const curve = new THREE.CatmullRomCurve3( Array.from(this.anchorPoints.children, (child) => child.position));
        
            const points = curve.getPoints( this.ARC_SEGMENTS );
            const geometry = new THREE.BufferGeometry().setFromPoints( points );
            this.curve.geometry = geometry;
        }

        this.anchorPoints.add(sphere);
        this.transformControls.attach(sphere);
    }

    removeCurve() {
        while (this.anchorPoints.children.length > 0) {
            this.anchorPoints.remove(this.anchorPoints.children[0]);
        }
        this.transformControls.detach();
        this.curve.geometry.dispose();
        this.curve.geometry = new THREE.BufferGeometry();

    }

}