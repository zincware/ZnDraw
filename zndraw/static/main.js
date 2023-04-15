import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';


// THREE.Cache.enabled = true;

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);

const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

const hemisphere_light = new THREE.HemisphereLight(0xffffff, 0x777777, 1);
scene.add(hemisphere_light);


const atoms = new THREE.Group();

const bonds = new THREE.Group();
const bonds_1 = new THREE.Group();
const bonds_2 = new THREE.Group();

let node1 = new THREE.Vector3();
let node2 = new THREE.Vector3();

// Helper Functions

function halfCylinderGeometry(pointX, pointY) {
	// Make the geometry (of "direction" length)
	var direction = new THREE.Vector3().subVectors(pointY, pointX);

	var geometry = new THREE.CylinderGeometry(0.15, 0.15, direction.length() / 2, 16);
	// // shift it so one end rests on the origin
	geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(0, direction.length() / 4, 0));
	// // rotate it the right way for lookAt to work
	geometry.applyMatrix4(new THREE.Matrix4().makeRotationX(THREE.MathUtils.degToRad(90)));

	return geometry;

}

function halfCylinderMesh(pointX, pointY, material) {
	// // Make a mesh with the geometry
	var geometry = halfCylinderGeometry(pointX, pointY);
	var mesh = new THREE.Mesh(geometry, material);
	// alignBetweenVectors(pointX, pointY, mesh);
	// // Position it where we want
	mesh.position.copy(pointX);
	// // And make it point to where we want
	mesh.lookAt(pointY);

	return mesh;
}

// Setup Scene

function addAtom(item) {
	const geometry = new THREE.SphereGeometry(item["radius"], 32, 16);
	const material = new THREE.MeshPhongMaterial({ color: item["color"] });
	const particle = new THREE.Mesh(geometry, material);
	atoms.add(particle);
	particle.position.set(...item["position"]);
	particle.callback = function () {
		let data = {
			"position": this.position,
		}
		fetch("atom/" + item["id"], {
			"method": "POST",
			"headers": { "Content-Type": "application/json" },
			"body": JSON.stringify(data),
		})
	}
}

function addBond(item) {
	atoms.children[item[0]].getWorldPosition(node1);
	atoms.children[item[1]].getWorldPosition(node2);

	const bond_1 = halfCylinderMesh(node1, node2, atoms.children[item[0]].material);
	const bond_2 = halfCylinderMesh(node2, node1, atoms.children[item[1]].material);

	bond_1.userData["id"] = item[0];
	bond_2.userData["id"] = item[1];

	bonds_1.add(bond_1);
	bonds_2.add(bond_2);

}

function cleanScene() {
	while (atoms.children.length > 0) {
		scene.remove(atoms.children.shift());
	};
	while (bonds_1.children.length > 0) {
		scene.remove(bonds_1.children.shift());
	};
	while (bonds_2.children.length > 0) {
		scene.remove(bonds_2.children.shift());
	};
}


function drawAtoms(obj) {
	cleanScene();

	obj["nodes"].forEach(function (item, index) {
		// console.log("Adding item " + index + " to scene(" + item + ")");
		addAtom(item);
	});

	scene.add(atoms);


	obj["edges"].forEach(function (item, index) {
		// console.log("Adding item " + index + " to scene(" + item + ")");
		addBond(item);

	});

	bonds.add(bonds_1);
	bonds.add(bonds_2);
	scene.add(bonds);

}

let obj = await (await fetch("xyz")).json();
drawAtoms(obj);

// interactions

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

let position = [];


async function onPointerDown(event) {

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
		console.log(intersects[i].object.position);
		if (position.length === 0) {
			position = await (await fetch("animation")).json();
		}
		cleanScene();
		drawAtoms(obj);

	}
}

// Async with onInterval can spawn many requests
async function onInterval(event) {
	// drawAtoms(obj);
	// drawAtoms();
	// let obj = await (await fetch("update")).json();
	// if (obj.hasOwnProperty('position')) {
	// 	// atoms.children[obj["id"]].position.set(...obj["position"])
	// 	//atoms.children[obj["id"]].material.color.set(0x00ff00);
	// 	cleanScene();
	// 	drawAtoms();

	// 	// const geometry = new THREE.SphereGeometry(obj["radius"], 32, 16);
	// 	// const material = new THREE.MeshPhongMaterial({ color: obj["color"] });
	// 	// const cube = new THREE.Mesh(geometry, material);
	// 	// atoms.add(cube);
	// 	// cube.position.set(...obj["position"]);
	// }
}

window.addEventListener('pointerdown', onPointerDown, false);
window.addEventListener('resize', onWindowResize, false);
// window.setInterval(onInterval, 1000);

function animate() {

	renderer.render(scene, camera);
	controls.update();

	if (position.length > 0) {
		let theRemovedElement = position.shift();
		theRemovedElement.forEach(function (item, index) {
			atoms.children[index].position.set(...item);
			// atoms.children[index].position.x += item[0];
			// atoms.children[index].position.y += item[1];
			// atoms.children[index].position.z += item[2];
		});
		console.log("Animation running")

		scene.updateMatrixWorld();

		for (let i = 0; i < bonds_1.children.length; i++) {
			// can't resize the cylinders
			
			let bond_1 = bonds_1.children[i];
			let bond_2 = bonds_2.children[i];
			
			atoms.children[bond_1.userData["id"]].getWorldPosition(node1);
			atoms.children[bond_2.userData["id"]].getWorldPosition(node2);

			let direction = new THREE.Vector3().subVectors(node1, node2);

			let scale = (direction.length() / 2) / bond_1.geometry.parameters.height;

			if (scale > 2) {
				bond_1.material.transparent = true ;
				bond_2.material.transparent = true ;
			}

			bond_1.scale.set(1, 1, scale);
			bond_2.scale.set(1, 1, scale);

			bond_1.position.copy(node1);
			bond_2.position.copy(node2);

			bond_1.lookAt(node2);
			bond_2.lookAt(node1);
		}
	}



	// animation loop

	requestAnimationFrame(animate);


}


animate();
