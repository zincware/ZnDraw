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

/**
 * ThreeJS variables
 */

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });

const hemisphereLight = new THREE.HemisphereLight(0xffffff, 0x777777, 0.1);
const spotLight = new THREE.SpotLight(0xffffff, 1, 0, Math.PI / 2);

const atomsGroup = new THREE.Group();

const bondsGroup = new THREE.Group();
const bondsGroup_1 = new THREE.Group();
const bondsGroup_2 = new THREE.Group();

let node1 = new THREE.Vector3();
let node2 = new THREE.Vector3();

const materials = {
	"MeshBasicMaterial": new THREE.MeshBasicMaterial({ color: "#ffa500" }),
	"MeshLambertMaterial": new THREE.MeshLambertMaterial({ color: "#ffa500" }),
	"MeshMatcapMaterial": new THREE.MeshMatcapMaterial({ color: "#ffa500" }),
	"MeshPhongMaterial": new THREE.MeshPhongMaterial({ color: "#ffa500" }),
	"MeshPhysicalMaterial": new THREE.MeshPhysicalMaterial({ color: "#ffa500" }),
	"MeshStandardMaterial": new THREE.MeshStandardMaterial({ color: "#ffa500" }),
	"MeshToonMaterial": new THREE.MeshToonMaterial({ color: "#ffa500" }),

};

/**
 * Three JS Setup
 */

renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);


spotLight.position.set(0, 0, 100);
scene.add(spotLight);

scene.add(hemisphereLight);

let config = {};


// some global variables
let frames = [];
let selected_ids = [];
let animation_frame = 0;
let scene_building = false;
let animation_running = true;
let data_loading = false;
let fps = [];

let keydown = { "shift": false, "ctrl": false, "alt": false, "c": false, "l": false };

/**
 * DOM variables
 */

const div_info = document.getElementById('info');
const div_loading = document.getElementById('loading');
const div_progressBar = document.getElementById('progressBar');
const div_bufferBar = document.getElementById('bufferBar');
const div_greyOut = document.getElementById('greyOut');
const div_lst_selected_ids = document.getElementById('lst_selected_ids');
const div_FPS = document.getElementById('FPS');
const div_n_particles = document.getElementById('n_particles');
const div_n_bonds = document.getElementById('n_bonds');
const div_help_container = document.getElementById('help_container');
const div_python_class_control = document.getElementById('python_class_control');

const o_selectAtoms = document.getElementById('selectAtoms');
const o_autoRestart = document.getElementById('autoRestart');
const o_animate = document.getElementById('animate');
const o_reset_selection = document.getElementById('reset_selection');
const o_hide_selection = document.getElementById('hide_selection');
const o_reset = document.getElementById('reset');
const o_max_fps = document.getElementById('max_fps');
const o_frames_per_post = document.getElementById('frames_per_post');
const o_sphere_plus = document.getElementById('sphere_plus');
const o_sphere_minus = document.getElementById('sphere_minus');
const o_bond_plus = document.getElementById('bond_plus');
const o_bond_minus = document.getElementById('bond_minus');
const o_resolution_plus = document.getElementById('resolution_plus');
const o_resolution_minus = document.getElementById('resolution_minus');
const o_materialSelect = document.getElementById('materialSelect');
const o_wireframe = document.getElementById('wireframe');
const o_spotLightIntensity = document.getElementById('spotLightIntensity');
const o_hemisphereLightIntensity = document.getElementById('hemisphereLightIntensity');
const o_help_btn = document.getElementById('help_btn');
const o_add_btn = document.getElementById('add_btn');
const o_newPythonClassBtn = document.getElementById('newPythonClassBtn');


// Helper Functions

async function load_config() {
	config = await (await fetch("config")).json();
	console.log(config)
}

async function update_materials() {
	for (const material in materials) {
		const option = document.createElement("option");
		option.text = material;
		option.value = material;
		o_materialSelect.appendChild(option);
	}
	o_materialSelect.value = "MeshPhongMaterial";
}

update_materials();

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

	const particle = new THREE.Mesh(geometry, materials[o_materialSelect.value].clone());
	atomsGroup.add(particle);
	particle.material.color.set(item["color"]);
	particle.material.wireframe = o_wireframe.checked;

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
	if (scene_building) {
		return;
	}
	const urls = ["atoms", "bonds"];
	animation_frame = step;
	console.log("Updating scene");

	div_loading.style.visibility = 'visible';
	div_greyOut.style.visibility = 'visible';

	// this is faster then doing it one by one
	const arrayOfResponses = await Promise.all(
		urls.map((url) =>
			fetch(url, {
				"method": "POST",
				"headers": { "Content-Type": "application/json" },
				"body": JSON.stringify(step)
			})
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

	if (intersects.length == 0) {
		return;
	}

	// for (let i = 0; i < intersects.length; i++) {
	let mesh = intersects[0].object;
	if (!keydown["shift"]) {
		selected_ids = [mesh.userData["id"]];
		mesh.material.color.set(0xffa500);
	} else {
		if (!selected_ids.includes(mesh.userData["id"])) {
			mesh.material.color.set(0xffa500);
			selected_ids.push(mesh.userData["id"]);
		} else {
			mesh.material.color.set(mesh.userData["color"]);
			selected_ids.splice(selected_ids.indexOf(mesh.userData["id"]), 1);
		}
	}
	// update colors here for better performance
	update_color_of_ids(selected_ids);
	await update_selection();
}

/**
 * We update the color of every atom in the scene
 * @param {list[int]} ids, the selected ids
 * @returns 
 */

async function update_color_of_ids(ids) {
	atomsGroup.children.forEach(function (mesh) {
		if (ids.includes(mesh.userData["id"])) {
			mesh.material.color.set(0xffa500);
		} else {
			mesh.material.color.set(mesh.userData["color"]);
		}
	});
	return ids;
}

async function update_selection() {
	console.log("Updating selection");
	div_lst_selected_ids.innerHTML = selected_ids.join(", ");
	await fetch("select", {
		"method": "POST",
		"headers": { "Content-Type": "application/json" },
		"body": JSON.stringify({ "selected_ids": selected_ids, "step": animation_frame }),
	}).then(response => response.json()).then(function (response_json) {
		if (response_json["updated"]) {
			update_color_of_ids(response_json["selected_ids"]);
			selected_ids = response_json["selected_ids"];
		}
		div_lst_selected_ids.innerHTML = response_json["selected_ids"].join(", ");
	});

}

async function getAnimationFrames() {
	if (data_loading) {
		console.log("Animation read already in progress");
		return;
	}
	console.log("Loading animation frames");

	data_loading = true;
	if (frames.length < 2) {
		fetch("load");
	}

	let step = frames.length;
	while (true) {

		let obj = await fetch("positions", {
			"method": "POST",
			"headers": { "Content-Type": "application/json" },
			"body": JSON.stringify({ "start": step, "stop": parseInt(o_frames_per_post.value) + step }),
		}).then(response => response.json()).then(function (response_json) {
			load_config();
			return response_json;
		});

		console.log("Read " + step + "-" + (parseInt(o_frames_per_post.value) + step) + " frames");
		if (Object.keys(obj).length === 0) {
			console.log("Animation read finished");
			break;
		}
		frames = frames.concat(obj);
		step += parseInt(o_frames_per_post.value);
	}
	data_loading = false;
}


if (config["animate"] === true) {
	div_info.innerHTML = "Reading file...";
	getAnimationFrames();
}
if (config["restart_animation"] === true) {
	o_autoRestart.checked = true;
}

/**
 * Event listeners
 */

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
	});
}

o_reset.onclick = function () {
	load_config();
	animation_frame = 0;
	build_scene(0);
	selected_ids = [];
	update_selection();
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

o_materialSelect.onchange = function () {
	build_scene(animation_frame);
}

o_wireframe.onchange = function () {
	build_scene(animation_frame);
}

o_spotLightIntensity.oninput = function () {
	document.getElementById("spotLightIntensity_output").value = o_spotLightIntensity.value;
	spotLight.intensity = o_spotLightIntensity.value;
}

o_hemisphereLightIntensity.oninput = function () {
	document.getElementById("hemisphereLightIntensity_output").value = o_hemisphereLightIntensity.value;
	hemisphereLight.intensity = o_hemisphereLightIntensity.value;
}

o_help_btn.onmouseover = function () {
	div_help_container.style.display = "block";
}

o_help_btn.onmouseout = function () {
	div_help_container.style.display = "none";
}

o_add_btn.onclick = function () {
	document.getElementById("add_class").style.display = "block";
}


/**
 * Helper function, move later
 * @param {} name 
 * @param {*} checked 
 * @returns 
 */
function createRadioElement(name, checked, id, properties) {
	var radioHtml = '<input class="form-check-input" type="radio" name="' + name + '"  id="' + id + '"';
	if (checked) {
		radioHtml += ' checked="checked"';
	}
	radioHtml += '/>';

	var radioFragment = document.createElement('div');
	radioFragment.classList.add("form-check");
	radioFragment.innerHTML = radioHtml;



	let function_container = document.createElement('div');
	function_container.classList.add("container-fluid", "bg-light", "rounded", "border", "border-primary");

	let function_container_label = document.createElement('h5');
	function_container_label.innerHTML = id;

	function_container.appendChild(function_container_label);

	let function_container_col = document.createElement('div');
	function_container_col.classList.add("row");

	let descriptions = document.createElement('div');
	descriptions.classList.add("col-sm-2");

	let values = document.createElement('div');
	values.classList.add("col-sm-1");

	let controllers = document.createElement('div');
	controllers.classList.add("col-sm-8");

	function_container_col.appendChild(descriptions);
	function_container_col.appendChild(values);
	function_container_col.appendChild(controllers);

	console.log(properties);

	Object.values(properties).forEach((item) => {
		console.log(item);
		let label = document.createElement('div');
		label.innerHTML = item["title"];
		let label_row = document.createElement('div');
		label_row.classList.add("row-sm");
		label_row.appendChild(label);

		let value = document.createElement('div');
		value.innerHTML = item["default"];
		let value_row = document.createElement('div');
		value_row.classList.add("row-sm");
		value_row.appendChild(value);

		descriptions.appendChild(label_row);
		values.appendChild(value_row);

		let controller = document.createElement('input');
		if (item["type"] == "integer") {
			controller.type = "range";
			controller.step = 1;
		} else if (item["type"] == "number") {
			controller.type = "range";
			controller.step = 0.1;
		} else if (item["type"] == "text") {
			controller.type = "text";
		} else {
			console.log("Unknown type: " + item["type"]);
		}
		controller.value = item["default"];
		controller.id = id + "_" + item["title"];

		if ("minimum" in item) {
			controller.min = item["minimum"];
		}
		if ("maximum" in item) {
			controller.max = item["maximum"];
		}

		let controller_row = document.createElement('div');
		controller_row.classList.add("row-sm");
		controller_row.appendChild(controller);

		controllers.appendChild(controller_row);

		controller.oninput = function () {
			value.innerHTML = controller.value;
		}

		controller.onchange = function () {
			// fetch with post 
			fetch("update_function_values", {
				"method": "POST",
				"headers": { "Content-Type": "application/json" },
				"body": JSON.stringify({
					"function_id": id,
					"property": item["title"],
					"value": controller.value
				})
			});
		}
	});
	function_container.appendChild(function_container_col);
	radioFragment.appendChild(function_container);

	return radioFragment;
}

o_newPythonClassBtn.onclick = function () {
	document.getElementById("add_class").style.display = "none";

	fetch("add_update_function", {
		"method": "POST",
		"headers": { "Content-Type": "application/json" },
		"body": JSON.stringify(document.getElementById("newPythonClass").value),
	}).then(response => response.json()).then(function (response_json) {
		// if not null alert
		if ("error" in response_json) {
			alert(response_json["error"]);
			stepError(response_json["error"]);
		} else {
			console.log(response_json);
			load_config();
		}
		return response_json;
	}).then(function (response_json) {
		console.log(response_json);
		div_python_class_control.appendChild(createRadioElement("flexRadioUpdateFunction", true, response_json["title"], response_json["properties"]));

		document.getElementById(response_json["title"]).onclick = function () {
			console.log("clicked");
			console.log(document.querySelector('input[name="flexRadioUpdateFunction"]:checked').id);
			fetch("/select_update_function/" + document.querySelector('input[name="flexRadioUpdateFunction"]:checked').id)
		};

	});


}

window.addEventListener("keydown", (event) => {
	if (event.isComposing || event.key === " ") {
		event.preventDefault();
		if (animation_running && animation_frame == frames.length - 1) {
			animation_frame = 0;
		} else {
			animation_running = !animation_running;
		}
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
		animation_frame = parseInt(Math.min(frames.length - 1, animation_frame + (frames.length / 10)));
	}
	if (event.isComposing || event.key === "ArrowDown") {
		animation_frame = parseInt(Math.max(0, animation_frame - (frames.length / 10)));
	}
	if (event.isComposing || event.key === "q") {
		getAnimationFrames();
	}
	if (event.isComposing || event.shiftKey) {
		keydown["shift"] = true;
	}
	if (event.isComposing || event.ctrlKey) {
		keydown["strg"] = true;
	}
	if (event.isComposing || event.altKey) {
		keydown["alt"] = true;
	}
	for (let key in keydown) {
		if (event.isComposing || event.key === key) {
			keydown[key] = true;
		}
	}
});

window.addEventListener("keyup", (event) => {
	if (event.isComposing || !event.shiftKey) {
		keydown["shift"] = false;
	}
	if (event.isComposing || !event.ctrlKey) {
		keydown["ctrl"] = false;
	}
	if (event.isComposing || !event.altKey) {
		keydown["alt"] = false;
	}
	for (let key in keydown) {
		if (event.isComposing || event.key === key) {
			keydown[key] = false;
		}
	}
});

window.addEventListener("keydown", (event) => {
	if (event.isComposing || event.key === "Enter") {
		div_info.innerHTML = "Processing...";

		fetch("update", {
			"method": "POST",
			"headers": { "Content-Type": "application/json" },
			"body": JSON.stringify({ "selected_ids": selected_ids, "step": animation_frame }),
		}).then((response) => getAnimationFrames());

		if (!data_loading) {
			getAnimationFrames();
		}
	}
});




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
		atomsGroup.getObjectByUserDataProperty("id", index).position.set(...item);
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

	fps.push(move_atoms_clock.getElapsedTime());

	if (fps.length > 10) {
		fps.shift();
	}
	if (!data_loading) {
		getAnimationFrames();
	}

	div_FPS.innerHTML = (1 / (fps.reduce((a, b) => a + b, 0) / fps.length)).toFixed(2);
	move_atoms_clock.start();
}

function centerCamera() {
	if (selected_ids.length === 0) {
		controls.target = new THREE.Vector3(0, 0, 0);
	} else {
		controls.target = atomsGroup.getObjectByUserDataProperty("id", selected_ids[0]).position.clone();
	}
}

function animate() {

	renderer.render(scene, camera);
	controls.update();

	if (frames.length > 0) {
		move_atoms();
	}
	if (keydown["c"]) {
		centerCamera();
	}
	if (keydown["l"]) {
		spotLight.position.x = camera.position.x;
		spotLight.position.y = camera.position.y;
		spotLight.position.z = camera.position.z;
	}

	// animation loop

	requestAnimationFrame(animate);
}

animate();
