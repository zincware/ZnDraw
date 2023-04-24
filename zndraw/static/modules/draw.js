import * as THREE from 'three';
import { TransformControls } from 'three/addons/controls/TransformControls.js';


let camera, renderer, scene, orbit;
let spline;
export const positions = [];
const anchorPoints = [];
export let control;

const ARC_SEGMENTS = 200;


export function init(_camera, _renderer, _scene, _orbit) {
    camera = _camera;
    renderer = _renderer;
    scene = _scene;
    orbit = _orbit;
    control = new TransformControls(camera, renderer.domElement);
    control.addEventListener('dragging-changed', function (event) {

        orbit.enabled = !event.value;

    });
    control.addEventListener('objectChange', function () {

        updateSplineOutline();

    });
    createCurve();
};

/**
 * Add an anchor point for drawing to the scene
 * @param {THREE.Vector3} position 
 * @returns 
 */
export function createAnchorPoint(position) {
    console.log("Creating anchor point at: ");
    console.log(position);
    let geometry = new THREE.SphereGeometry(0.2, 32, 32);
    let material = new THREE.MeshBasicMaterial({ color: "#000000" });
    let sphere = new THREE.Mesh(geometry, material);
    sphere.position.copy(position);

    if (positions.length > 0) {
        // remove last entry from positions and remplace with clone from last_position
        if (positions.length === 1) { // we are currently adding the second point
            scene.add(spline.mesh);
        }
        let last_position = positions.pop();
        positions.push(last_position.clone());
    } else {
        scene.add(control);
    }
    positions.push(sphere.position);
    // control detach first?
    control.attach(sphere);
    scene.add(sphere);
    anchorPoints.push(sphere);

    updateSplineOutline();
    renderer.render(scene, camera);
}

function createCurve() {

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(ARC_SEGMENTS * 3), 3));

    let curve = new THREE.CatmullRomCurve3(positions);
    curve.curveType = 'catmullrom';
    curve.mesh = new THREE.Line(geometry.clone(), new THREE.LineBasicMaterial({
        color: 0xff0000,
        opacity: 0.35
    }));
    curve.mesh.castShadow = true;
    spline = curve;
}

function updateSplineOutline() {

    const point = new THREE.Vector3();

    const splineMesh = spline.mesh;
    const position = splineMesh.geometry.attributes.position;

    for (let i = 0; i < ARC_SEGMENTS; i++) {

        const t = i / (ARC_SEGMENTS - 1);
        spline.getPoint(t, point);
        position.setXYZ(i, point.x, point.y, point.z);

    }

    position.needsUpdate = true;
    renderer.render(scene, camera);

}

export function reset() {
    for (let i = 0; i < anchorPoints.length; i++) {
        scene.remove(anchorPoints[i]);
    }
    positions.length = 0;
    scene.remove(spline.mesh);
    scene.remove(control);
    updateSplineOutline();
}