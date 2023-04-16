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
let animation_running = true;
const div_info = document.getElementById('info');
const div_loading = document.getElementById('loading');
const div_progressBar = document.getElementById('progressBar');
const div_bufferBar = document.getElementById('bufferBar');
const div_greyOut = document.getElementById('greyOut');
const div_lst_selected_ids = document.getElementById('lst_selected_ids');
const div_FPS = document.getElementById('FPS');
const div_n_particles = document.getElementById('n_particles');
const div_n_bonds = document.getElementById('n_bonds');

const o_selectAtoms = document.getElementById('selectAtoms');
const o_autoRestart = document.getElementById('autoRestart');
const o_animate = document.getElementById('animate');
const o_reset_selection = document.getElementById('reset_selection');
const o_hide_selection = document.getElementById('hide_selection');
const o_reset = document.getElementById('reset');
const o_max_fps = document.getElementById('max_fps');
const o_frame_buffer = document.getElementById('frame_buffer');
const o_frames_per_post = document.getElementById('frames_per_post');
const o_sphere_plus = document.getElementById('sphere_plus');
const o_sphere_minus = document.getElementById('sphere_minus');
const o_bond_plus = document.getElementById('bond_plus');
const o_bond_minus = document.getElementById('bond_minus');
const o_resolution_plus = document.getElementById('resolution_plus');
const o_resolution_minus = document.getElementById('resolution_minus');


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
	particle.userData["bond_ids"] = [];
}

function addBond(item) {
	atomsGroup.children[item[0]].getWorldPosition(node1);
	atomsGroup.children[item[1]].getWorldPosition(node2);

	const bond_1 = halfCylinderMesh(node1, node2, atomsGroup.children[item[0]].material);
	const bond_2 = halfCylinderMesh(node2, node1, atomsGroup.children[item[1]].material);

	bond_1.userData["atom_id"] = item[0];
	bond_2.userData["atom_id"] = item[1];

	atomsGroup.children[item[0]].userData["bond_ids"].push(bond_1.id);
	atomsGroup.children[item[0]].userData["bond_ids"].push(bond_2.id);
	atomsGroup.children[item[1]].userData["bond_ids"].push(bond_1.id);
	atomsGroup.children[item[1]].userData["bond_ids"].push(bond_2.id);


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
	animation_frame = step;
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
	await update_selection();
	scene_building = false;

	div_n_particles.innerHTML = atomsGroup.children.length;
	div_n_bonds.innerHTML = bondsGroup_1.children.length;

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
	await update_selection();
}

async function update_selection() {
	fetch("select", {
		"method": "POST",
		"headers": { "Content-Type": "application/json" },
		"body": JSON.stringify(selected_ids),
	})

	div_lst_selected_ids.innerHTML = selected_ids.join(", ");
}

async function getAnimationFrames() {
	frames = [];

	await fetch("atoms/1")
	load_config();

	let step = 0;
	while (true) {
		let obj = await (await fetch("atoms/" + step + "&" + (parseInt(o_frames_per_post.value) + step))).json();
		if (Object.keys(obj).length === 0) {
			console.log("Animation read finished");
			break;
		}
		frames = frames.concat(obj);  // TODO: handle multiple frames at once
		step += parseInt(o_frames_per_post.value);
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

o_reset_selection.onclick = function () {
	selected_ids.forEach(function (item, index) {
		let mesh = atomsGroup.getObjectByUserDataProperty("id", item);
		mesh.material.color.set(mesh.userData["color"]);
	});
	selected_ids = [];
	update_selection();
}

o_hide_selection.onclick = function () {
	selected_ids.forEach(function (item, index) {
		let mesh = atomsGroup.getObjectByUserDataProperty("id", item);
		mesh.visible = false;

		mesh.userData["bond_ids"].forEach(function (item, index) {
			bondsGroup_1.getObjectById(item).visible = false;
			bondsGroup_2.getObjectById(item).visible = false;
		});

		// bondsGroup_1.getObjectById(item).visible = false;
		// bondsGroup_2.getObjectByUserDataProperty("id", item).visible = false;
	});
}

o_reset.onclick = function () {
	load_config();
	build_scene(0);
	selected_ids = [];
	update_selection();
	camera.position.z = 50;
}

o_sphere_plus.onclick = function () {
	config["sphere_size"] += 0.1;
	build_scene(animation_frame);
}

o_sphere_minus.onclick = function () {
	config["sphere_size"] -= 0.1;
	build_scene(animation_frame);
}

o_bond_plus.onclick = function () {
	config["bond_size"] += 0.1;
	build_scene(animation_frame);
}

o_bond_minus.onclick = function () {
	config["bond_size"] -= 0.1;
	build_scene(animation_frame);
}

o_resolution_plus.onclick = function () {
	config["resolution"] += 1;
	build_scene(animation_frame);
}

o_resolution_minus.onclick = function () {
	config["resolution"] -= 1;
	build_scene(animation_frame);
}

window.addEventListener("keydown", (event) => {
	if (event.isComposing || event.key === " ") {
		event.preventDefault();
		animation_running = !animation_running;
	}
	if (event.isComposing || event.key === "ArrowLeft") {
		event.preventDefault();
		animation_frame = Math.max(0, animation_frame - 1);
	}
	if (event.isComposing || event.key === "ArrowRight") {
		event.preventDefault();
		animation_frame = Math.min(frames.length - 1, animation_frame + 1);
	}
	if (event.isComposing || event.key === "ArrowUp") {
		animation_frame = frames.length - 1;
	}
	if (event.isComposing || event.key === "ArrowDown") {
		animation_frame = 0;
	}
});

if (config["update_function"] !== null) {
	window.addEventListener("keydown", (event) => {
		if (event.isComposing || event.key === "Enter") {
			div_info.innerHTML = "Processing...";
			fetch("update/" + animation_frame).then((response) => getAnimationFrames());
		}
	});
}



let move_atoms_clock = new THREE.Clock();


function move_atoms() {
	if (frames.length === 0) {
		return;
	}
	if (move_atoms_clock.getElapsedTime() < (1 / o_max_fps.value)) {
		return;
	}
	if (scene_building === true) {
		return;
	}
	console.log("Animation (" + animation_frame + "/" + (frames.length - 1) + ")");
	if (frames.length < parseInt(o_frame_buffer.value)) {
		div_info.innerHTML = "Buffering...";
		div_progressBar.style.width = (((frames.length - 1) / parseInt(o_frame_buffer.value)) * 100).toFixed(2) + "%";
		return;
	}
	div_progressBar.style.visibility = "hidden";

	if (animation_running === true) {
		if (animation_frame < frames.length - 1) {
			animation_frame += 1;
		} else if (o_autoRestart.checked === true) {
			animation_frame = 0;
		}
	}
	if (frames.length < animation_frame) {
		// waiting for async call to finish
		return;
	}
	if (config["total_frames"] > 0) {
		div_progressBar.style.visibility = "visible";
		div_progressBar.style.width = ((animation_frame / config["total_frames"]) * 100).toFixed(2) + "%";
		div_bufferBar.style.width = (((frames.length - 1) / config["total_frames"]) * 100).toFixed(2) + "%";
	}

	if (frames[animation_frame].length != atomsGroup.children.length) {
		// we need to update the scene
		build_scene(animation_frame);
		scene_building = true;
		return; // we need to wait for the scene to be updated
	}

	frames[animation_frame].forEach(function (item, index) {
		atomsGroup.getObjectByUserDataProperty("id", item["id"]).position.set(...item["position"]);
		// atomsGroup.children[item["id"]].position.set(...item["position"]);
	});

	div_info.innerHTML = "Frame (" + animation_frame + "/" + (frames.length - 1) + ")";

	if (config["bond_size"] > 0) {
		scene.updateMatrixWorld();

		for (let i = 0; i < bondsGroup_1.children.length; i++) {
			let bond_1 = bondsGroup_1.children[i];
			let bond_2 = bondsGroup_2.children[i];

			atomsGroup.getObjectByUserDataProperty("id", bond_1.userData["atom_id"]).getWorldPosition(node1);
			atomsGroup.getObjectByUserDataProperty("id", bond_2.userData["atom_id"]).getWorldPosition(node2);

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
	let fps = 1 / move_atoms_clock.getElapsedTime();
	div_FPS.innerHTML = fps.toFixed(2) + " FPS";
	move_atoms_clock.start();
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
