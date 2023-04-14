import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';


// THREE.Cache.enabled = true;

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);

const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

let obj = await (await fetch("xyz")).json();

obj["nodes"].forEach(function (item, index) {
	console.log("Adding item " + index + " to scene(" + item + ")");

	const geometry = new THREE.SphereGeometry(item["radius"], 32, 16);
	const material = new THREE.MeshBasicMaterial({ color: item["color"] });
	const cube = new THREE.Mesh(geometry, material);
	scene.add(cube);
	cube.position.set(item["x"], item["y"], item["z"]);
});

obj["edges"].forEach(function (item, index) {
	console.log("Adding item " + index + " to scene(" + item + ")");

	// const geometry = new THREE.CylinderGeometry(obj["radius"], obj["radius"], item["height"], 32);
	const geometry = new THREE.CylinderGeometry(0.1, 0.1, 1.0, 32);

	// const geometry = new THREE.SphereGeometry(item["radius"], 32, 16);
	const material = new THREE.MeshBasicMaterial({ color: "white" });
	const cylinder = new THREE.Mesh(geometry, material);
	scene.add(cylinder);
	cylinder.position.set(...item["position"]); // spread operator

	const quaternion = new THREE.Quaternion()
	const cylinderUpAxis = new THREE.Vector3( 0, 1, 0 )
	const stickAxis = new THREE.Vector3(...item["axis"])
	quaternion.setFromUnitVectors(cylinderUpAxis, stickAxis)
	cylinder.applyQuaternion(quaternion)
});


camera.position.z = 50;

const controls = new OrbitControls(camera, renderer.domElement);

// camera.position.set( 0, 20, 100 );
controls.update();

function animate() {

	requestAnimationFrame(animate);

	// required if controls.enableDamping or controls.autoRotate are set to true
	controls.update();

	renderer.render(scene, camera);

}

// function animate() {
// 	requestAnimationFrame( animate );

// 	cube.rotation.x += 0.01;
// 	cube.rotation.y += 0.01;

// 	renderer.render( scene, camera );
// }

animate();