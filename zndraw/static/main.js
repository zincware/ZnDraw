import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';


THREE.Object3D.prototype.getObjectByUserDataProperty = function (name, value) {

	if (this.userData[name] === value) return this;

	for (var i = 0, l = this.children.length; i < l; i++) {

		var child = this.children[i];
		var object = child.getObjectByUserDataProperty(name, value);

		if (object !== undefined) {

			return object;

		}

	}

	return undefined;

}

// THREE.Cache.enabled = true;

let config = {};

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

// some global variables
let frames = [];
let selected_ids = [];
let animation_frame = 0;
let scene_building = false;
const div_info = document.getElementById('info');
const div_loading = document.getElementById('loading');
const div_progressBar = document.getElementById('progressBar');
const div_bufferBar = document.getElementById('bufferBar');
const div_greyOut = document.getElementById('greyOut');
const o_selectAtoms = document.getElementById('selectAtoms');
const o_autoRestart = document.getElementById('autoRestart');
const o_animate = document.getElementById('animate');


// Helper Functions

async function load_config() {
	config = await (await fetch("config")).json();
	console.log(config)
}

function halfCylinderGeometry(pointX, pointY) {
	// Make the geometry (of "direction" length)
	var direction = new THREE.Vector3().subVectors(pointY, pointX);

	var geometry = new THREE.CylinderGeometry(0.15 * config["bond_size"], 0.15 * config["bond_size"], direction.length() / 2, config["resolution"] * 2);
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
	const geometry = new THREE.SphereGeometry(item["radius"] * config["sphere_size"], config["resolution"] * 4, config["resolution"] * 2);
	const material = new THREE.MeshPhongMaterial({ color: item["color"] });
	const particle = new THREE.Mesh(geometry, material);
	atomsGroup.add(particle);
	particle.position.set(...item["position"]);
	particle.userData["id"] = item["id"];
	particle.userData["color"] = item["color"];
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
	if (config["bond_size"] > 0) {

		bonds.forEach(function (item, index) {
			// console.log("Adding item " + index + " to scene(" + item + ")");
			addBond(item);

		});

		bondsGroup.add(bondsGroup_1);
		bondsGroup.add(bondsGroup_2);
		scene.add(bondsGroup);

	}


}

async function build_scene(step) {
	const urls = ["atoms/" + step, "bonds/" + step];
	console.log("Updating scene");

	div_loading.style.visibility = 'visible';
	div_greyOut.style.visibility = 'visible';

	// this is faster then doing it one by one
	const arrayOfResponses = await Promise.all(
		urls.map((url) =>
			fetch(url)
				.then((res) => res.json())
		)
	);

	drawAtoms(arrayOfResponses[0], arrayOfResponses[1]);
	selected_ids = [];
	scene_building = false;

	div_greyOut.style.visibility = 'hidden';
	div_loading.style.visibility = 'hidden';
}
await load_config();

build_scene(0);

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

		// mesh.callback();

		if (selected_ids.includes(mesh.userData["id"])) {
			mesh.material.color.set(mesh.userData["color"]);
			let index = selected_ids.indexOf(mesh.userData["id"]);
			if (index !== -1) {
				selected_ids.splice(index, 1);
			}
		} else {
			intersects[i].object.material.color.set(0xffa500);
			selected_ids.push(intersects[i].object.userData["id"]);
		};
	}

	fetch("select", {
		"method": "POST",
		"headers": { "Content-Type": "application/json" },
		"body": JSON.stringify(selected_ids),
	})
}

async function getAnimationFrames() {
	frames = [];

	await fetch("atoms/1")
	load_config();

	let step = 0;
	while (true) {
		let obj = await (await fetch("atoms/" + step + "&" + (step + config["frames_per_post"]))).json();
		if (Object.keys(obj).length === 0) {
			console.log("Animation read finished");
			break;
		}
		frames = frames.concat(obj);  // TODO: handle multiple frames at once
		step += config["frames_per_post"];
	}
}


if (config["animate"] === true) {
	div_info.innerHTML = "Reading file...";
	getAnimationFrames();
}
if (config["restart_animation"] === true) {
	o_autoRestart.checked = true;
}

console.log(config);

window.addEventListener('pointerdown', onPointerDown, false);
window.addEventListener('resize', onWindowResize, false);

o_selectAtoms.onclick = function () {
	if (o_selectAtoms.checked) {
		console.log("Selecting atoms");
		window.addEventListener('pointerdown', onPointerDown, false);
	}
	else {
		console.log("Deselecting atoms");
		window.removeEventListener('pointerdown', onPointerDown, false);
	}
}
o_animate.onclick = function () {
	div_info.innerHTML = "Reading file...";
	getAnimationFrames();
}


window.addEventListener("keydown", (event) => {
	if (event.isComposing || event.key === " ") {
		animation_frame = 0;
	}
});

if (config["update_function"] !== null) {
	window.addEventListener("keydown", (event) => {
		if (event.isComposing || event.key === "Enter") {
			fetch("update/" + animation_frame).then((response) => getAnimationFrames());
		}
	});
}



let clock = new THREE.Clock();


function move_atoms() {
	if (frames.length === 0) {
		return;
	}
	if (clock.getElapsedTime() < (1 / config["max_fps"])) {
		return;
	}
	console.log("Animation (" + animation_frame + "/" + frames.length + ")");
	if (frames.length < config["frame_buffer"]) {
		div_info.innerHTML = "Buffering...";
		div_progressBar.style.width = ((frames.length / config["frame_buffer"]) * 100).toFixed(2) + "%";
		return;
	}
	div_progressBar.style.visibility = "hidden";

	if (animation_frame < frames.length - 1) {
		animation_frame += 1;
	} else if (o_autoRestart.checked === true) {
		animation_frame = 0;
	}
	if (frames.length < animation_frame) {
		// waiting for async call to finish
		return;
	}
	if (config["total_frames"] > 0) {
		div_progressBar.style.visibility = "visible";
		div_progressBar.style.width = ((animation_frame / config["total_frames"]) * 100).toFixed(2) + "%";
		div_bufferBar.style.width = ((frames.length / config["total_frames"]) * 100).toFixed(2) + "%";
	}

	if (frames[animation_frame].length != atomsGroup.children.length) {
		// we need to update the scene
		if (scene_building === true) {
			return;
		}
		build_scene(animation_frame);
		scene_building = true;
		return; // we need to wait for the scene to be updated
	}

	frames[animation_frame].forEach(function (item, index) {
		atomsGroup.getObjectByUserDataProperty("id", item["id"]).position.set(...item["position"]);
		// atomsGroup.children[item["id"]].position.set(...item["position"]);
	});

	div_info.innerHTML = "Frame (" + animation_frame + "/" + frames.length + ")";

	if (config["bond_size"] > 0) {
		scene.updateMatrixWorld();

		for (let i = 0; i < bondsGroup_1.children.length; i++) {
			// can't resize the cylinders

			let bond_1 = bondsGroup_1.children[i];
			let bond_2 = bondsGroup_2.children[i];

			atomsGroup.getObjectByUserDataProperty("id", bond_1.userData["id"]).getWorldPosition(node1);
			atomsGroup.getObjectByUserDataProperty("id", bond_2.userData["id"]).getWorldPosition(node2);

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
	clock.start();
}

function animate() {

	renderer.render(scene, camera);
	controls.update();

	if (frames.length > 0) {
		move_atoms();
	}

	// animation loop

	requestAnimationFrame(animate);


}


animate();
