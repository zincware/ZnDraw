import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';


// THREE.Cache.enabled = true;

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);

const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

let obj = await (await fetch("xyz")).json();

const atoms = new THREE.Group();

obj["nodes"].forEach(function (item, index) {
	// console.log("Adding item " + index + " to scene(" + item + ")");

	const geometry = new THREE.SphereGeometry(item["radius"], 32, 16);
	const material = new THREE.MeshBasicMaterial({ color: item["color"] });
	const cube = new THREE.Mesh(geometry, material);
	atoms.add(cube);
	cube.position.set(item["x"], item["y"], item["z"]);
	cube.callback = function () {
		let data = {
			"atom": item,
		}
		fetch("atom/" + index, {
			"method": "POST",
			"headers": { "Content-Type": "application/json" },
			"body": JSON.stringify(data),
		})
	}
	cube.radius += 0.1;

});

scene.add(atoms);

const bonds = new THREE.Group();

obj["edges"].forEach(function (item, index) {
	// console.log("Adding item " + index + " to scene(" + item + ")");

	// const geometry = new THREE.CylinderGeometry(obj["radius"], obj["radius"], item["height"], 32);
	const geometry = new THREE.CylinderGeometry(0.1, 0.1, 1.0, 32);

	// const geometry = new THREE.SphereGeometry(item["radius"], 32, 16);
	const material = new THREE.MeshBasicMaterial({ color: "white" });
	const cylinder = new THREE.Mesh(geometry, material);
	bonds.add(cylinder);
	cylinder.position.set(...item["position"]); // spread operator

	const quaternion = new THREE.Quaternion()
	const cylinderUpAxis = new THREE.Vector3(0, 1, 0)
	const stickAxis = new THREE.Vector3(...item["axis"])
	quaternion.setFromUnitVectors(cylinderUpAxis, stickAxis)
	cylinder.applyQuaternion(quaternion)
});

scene.add(bonds);


camera.position.z = 50;

const controls = new OrbitControls(camera, renderer.domElement);

// camera.position.set( 0, 20, 100 );
controls.update();

function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight
    camera.updateProjectionMatrix()
    renderer.setSize(window.innerWidth, window.innerHeight)
    renderer.render(scene, camera);
}

const raycaster = new THREE.Raycaster();
const pointer = new THREE.Vector2();

function onPointerDown(event) {

	// calculate pointer position in normalized device coordinates
	// (-1 to +1) for both components
	// event.preventDefault(); # this doesn't work

	pointer.x = (event.clientX / window.innerWidth) * 2 - 1;
	pointer.y = - (event.clientY / window.innerHeight) * 2 + 1;

	// update the picking ray with the camera and pointer position
	raycaster.setFromCamera(pointer, camera);

	// calculate objects intersecting the picking ray
	const intersects = raycaster.intersectObjects(atoms.children);

	for (let i = 0; i < intersects.length; i++) {

		intersects[i].object.material.color.set(0xff0000);
		intersects[i].object.callback();

	}
}

window.addEventListener('pointerdown', onPointerDown, false);
window.addEventListener('resize', onWindowResize, false)


function animate() {

	requestAnimationFrame(animate);
	renderer.render(scene, camera);
	controls.update();

}


animate();