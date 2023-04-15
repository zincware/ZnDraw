import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';


// THREE.Cache.enabled = true;

const config = await (await fetch("config")).json();

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);

const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

const hemisphere_light = new THREE.HemisphereLight(0xffffff, 0x777777, 1);
scene.add(hemisphere_light);


const atomsGroup = new THREE.Group();

const bondsGroup = new THREE.Group();
const bondsGroup_1 = new THREE.Group();
const bondsGroup_2 = new THREE.Group();

let node1 = new THREE.Vector3();
let node2 = new THREE.Vector3();

// Helper Functions

function halfCylinderGeometry(pointX, pointY) {
	// Make the geometry (of "direction" length)
	var direction = new THREE.Vector3().subVectors(pointY, pointX);

	var geometry = new THREE.CylinderGeometry(0.15 * config["bond_size"], 0.15 * config["bond_size"], direction.length() / 2, 16);
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
	const geometry = new THREE.SphereGeometry(item["radius"] * config["sphere_size"], 32, 16);
	const material = new THREE.MeshPhongMaterial({ color: item["color"] });
	const particle = new THREE.Mesh(geometry, material);
	atomsGroup.add(particle);
	particle.position.set(...item["position"]);
	particle.userData["id"] = item["id"];
	particle.userData["color"] = item["color"];
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
	atomsGroup.children[item[0]].getWorldPosition(node1);
	atomsGroup.children[item[1]].getWorldPosition(node2);

	const bond_1 = halfCylinderMesh(node1, node2, atomsGroup.children[item[0]].material);
	const bond_2 = halfCylinderMesh(node2, node1, atomsGroup.children[item[1]].material);

	bond_1.userData["id"] = item[0];
	bond_2.userData["id"] = item[1];

	bondsGroup_1.add(bond_1);
	bondsGroup_2.add(bond_2);

}

function cleanScene() {
	while (atomsGroup.children.length > 0) {
		scene.remove(atomsGroup.children.shift());
	};
	while (bondsGroup_1.children.length > 0) {
		scene.remove(bondsGroup_1.children.shift());
	};
	while (bondsGroup_2.children.length > 0) {
		scene.remove(bondsGroup_2.children.shift());
	};
}


function drawAtoms(atoms, bonds) {
	cleanScene();

	atoms.forEach(function (item, index) {
		// console.log("Adding item " + index + " to scene(" + item + ")");
		addAtom(item);
	});

	scene.add(atomsGroup);


	bonds.forEach(function (item, index) {
		// console.log("Adding item " + index + " to scene(" + item + ")");
		addBond(item);

	});

	bondsGroup.add(bondsGroup_1);
	bondsGroup.add(bondsGroup_2);
	scene.add(bondsGroup);

}

let atoms = await (await fetch("atoms")).json();
let bonds = await (await fetch("bonds")).json();
drawAtoms(atoms, bonds);

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

let selected_ids = [];

async function onPointerDown(event) {

	// calculate pointer position in normalized device coordinates
	// (-1 to +1) for both components
	// event.preventDefault(); # this doesn't work

	pointer.x = (event.clientX / window.innerWidth) * 2 - 1;
	pointer.y = - (event.clientY / window.innerHeight) * 2 + 1;

	// update the picking ray with the camera and pointer position
	raycaster.setFromCamera(pointer, camera);

	// calculate objects intersecting the picking ray
	const intersects = raycaster.intersectObjects(atomsGroup.children);

	for (let i = 0; i < intersects.length; i++) {

		let mesh = intersects[i].object;

		mesh.callback();

		if (selected_ids.includes(mesh.userData["id"])) {
			mesh.material.color.set(mesh.userData["color"]);
			let index = selected_ids.indexOf(mesh.userData["id"]);
			if (index !== -1) {
				selected_ids.splice(index, 1);
			}
		} else {

			intersects[i].object.material.color.set(0xff0000);
			console.log(intersects[i].object.position);
			selected_ids.push(intersects[i].object.userData["id"]);
		};

		// if (position.length === 0) {
		// 	position = await (await fetch("animation")).json();
		// }
		// cleanScene();
		// drawAtoms(obj);

	}
}

async function getAnimationFrames(url) {
	position = [];
	while (true) {
		let obj = await (await fetch(url)).json();
		if (Object.keys(obj).length === 0) {
			console.log("Animation read finished");
			break;
		}
		position = position.concat(obj);
	}
}

console.log(config);

if (config["animate"] === true) {
	getAnimationFrames("animation");
}

window.addEventListener('pointerdown', onPointerDown, false);
window.addEventListener('resize', onWindowResize, false);

window.addEventListener("keydown", (event) => {
	if (event.isComposing || event.keyCode === 229) {
		return;
	}
	if (config["update_function"] !== "") {
		getAnimationFrames("update");
	}
});

let animation_frame = 0;
let clock = new THREE.Clock();

function move_atoms() {
	if (animation_frame < position.length - 1) {
		animation_frame += 1;
	} else {
		animation_frame = 0;
	}
	if (clock.getElapsedTime() < 1 / config["max_fps"]) {
		return;
	}
	clock.start();

	position[animation_frame].forEach(function (item, index) {
		atomsGroup.children[index].position.set(...item);
	});
	console.log("Animation running")

	scene.updateMatrixWorld();

	for (let i = 0; i < bondsGroup_1.children.length; i++) {
		// can't resize the cylinders

		let bond_1 = bondsGroup_1.children[i];
		let bond_2 = bondsGroup_2.children[i];

		atomsGroup.children[bond_1.userData["id"]].getWorldPosition(node1);
		atomsGroup.children[bond_2.userData["id"]].getWorldPosition(node2);

		let direction = new THREE.Vector3().subVectors(node1, node2);

		let scale = (direction.length() / 2) / bond_1.geometry.parameters.height;

		if (scale > 1.5) {
			scale = 0.0;
		}

		bond_1.scale.set(1, 1, scale);
		bond_2.scale.set(1, 1, scale);

		bond_1.position.copy(node1);
		bond_2.position.copy(node2);

		bond_1.lookAt(node2);
		bond_2.lookAt(node1);
	}
}

function animate() {

	renderer.render(scene, camera);
	controls.update();

	if (position.length > 0) {
		move_atoms();
	}

	// animation loop

	requestAnimationFrame(animate);


}


animate();
