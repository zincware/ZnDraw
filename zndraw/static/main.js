import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { particleGroup, materials, drawAtoms, speciesMaterial, countBonds, getAtomById, updateParticlePositions, printIndices } from './modules/particles.js';
import { keydown } from './modules/keypress.js';

// THREE.Cache.enabled = true;

/**
 * ThreeJS variables
 */

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });

const hemisphereLight = new THREE.HemisphereLight(0xffffff, 0x777777, 0.1);
const spotLight = new THREE.SpotLight(0xffffff, 1, 0, Math.PI / 2);



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
let displayed_frame = 0;
let scene_building = false;
let animation_running = true;
let data_loading = false;
let fps = [];


/**
 * DOM variables
 */

const div_info = document.getElementById('info');
const div_loading = document.getElementById('atom-spinner');
const div_progressBar = document.getElementById('progressBar');
const div_bufferBar = document.getElementById('bufferBar');
const div_greyOut = document.getElementById('greyOut');
const div_lst_selected_ids = document.getElementById('lst_selected_ids');
const div_FPS = document.getElementById('FPS');
const div_n_particles = document.getElementById('n_particles');
const div_n_bonds = document.getElementById('n_bonds');
const div_help_container = document.getElementById('help_container');
const div_python_class_control = document.getElementById('python_class_control');

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

const addModifierModal = new bootstrap.Modal(document.getElementById("addModifierModal"));
const addSceneModifier = document.getElementById("addSceneModifier");


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



// Setup Scene

let build_scene_cache = {};

async function build_scene(step) {
	if (scene_building) {
		return;
	}
	console.log("Updating scene");

	div_loading.style.visibility = 'visible';
	div_greyOut.style.visibility = 'visible';

	animation_frame = step;

	let arrayOfResponses = [];
	if (build_scene_cache.hasOwnProperty(step)) {
		console.log("Using cached scene");
		arrayOfResponses = build_scene_cache[step];
	} else {
		// this is faster then doing it one by one

		arrayOfResponses = await (await fetch("graph", {
			"method": "POST",
			"headers": { "Content-Type": "application/json" },
			"body": JSON.stringify(step)
		})).json();
		console.log(arrayOfResponses);
		build_scene_cache[step] = arrayOfResponses;
	}

	drawAtoms(arrayOfResponses["nodes"], arrayOfResponses["edges"], config, scene);
	selected_ids = [];
	await update_selection();
	scene_building = false;

	div_n_particles.innerHTML = particleGroup.children.length;

	div_n_bonds.innerHTML = countBonds();

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
	const intersects = raycaster.intersectObjects(particleGroup.children, true);

	if (intersects.length == 0) {
		return;
	}

	// for (let i = 0; i < intersects.length; i++) {
	let mesh = intersects[0].object;
	if (!keydown["shift"]) {
		if (selected_ids.includes(mesh.userData["id"])) {
			selected_ids = [];
		} else {
			selected_ids = [mesh.userData["id"]];
		}
	} else {
		if (!selected_ids.includes(mesh.userData["id"])) {
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
	particleGroup.children.forEach(function (atom_grp) {
		let mesh = atom_grp.children[0];
		if (ids.includes(mesh.userData["id"])) {
			let material = speciesMaterial(document.getElementById('materialSelect').value, 0xffa500, document.getElementById('wireframe').checked);
			atom_grp.children.forEach(function (item) {
				item.material = material;
			});
		} else {
			let material = speciesMaterial(document.getElementById('materialSelect').value, mesh.userData["color"], document.getElementById('wireframe').checked);
			atom_grp.children.forEach(function (item) {
				item.material = material;
			});
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
		"body": JSON.stringify({ "selected_ids": selected_ids, "step": animation_frame, "method": document.getElementById("selection-method").value }),
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

document.getElementById("selection-method").onclick = function () {
	if (document.getElementById("selection-method").value == "none") {
		console.log("Deselecting atoms");
		window.removeEventListener('pointerdown', onPointerDown, false);
	} else {
		console.log("Selecting atoms");
		window.addEventListener('pointerdown', onPointerDown, false);
	}

}

o_animate.onclick = function () {
	div_info.innerHTML = "Reading file...";
	getAnimationFrames();
}


o_reset_selection.onclick = function () {
	selected_ids = [];
	update_color_of_ids(selected_ids);
	update_selection();
}


// o_hide_selection.onclick = function () {
// 	selected_ids.forEach(function (item, index) {
// 		let mesh = particleGroup.getObjectByUserDataProperty("id", item);
// 		mesh.visible = false;

// 		mesh.userData["bond_ids"].forEach(function (item, index) {
// 			bondsXGroup.children[0].getObjectById(item).visible = false;
// 			bondsXGroup.children[1].getObjectById(item).visible = false;
// 		});
// 	});
// }

o_reset.onclick = function () {
	load_config();
	animation_frame = 0;
	build_scene(0);
	selected_ids = [];
	update_selection();
}

document.getElementById("sphereRadius").onchange = function () {
	let radius = parseFloat(document.getElementById("sphereRadius").value);
	let particleGeometry = particleGroup.children[0].children[0].geometry;
	let scale = radius / particleGeometry.boundingSphere.radius;
	particleGeometry.scale(scale, scale, scale);
};

document.getElementById("sphereRadius").oninput = function () {
	let radius = parseFloat(document.getElementById("sphereRadius").value);
	document.getElementById("sphereRadiusLabel").innerHTML = "Sphere radius: " + radius;
};

document.getElementById("bondDiameter").oninput = function () {
	let radius = parseFloat(document.getElementById("bondDiameter").value);
	document.getElementById("bondDiameterLabel").innerHTML = "Bond diameter: " + radius;
};

document.getElementById("bondDiameter").onchange = function () {
	let diameter = parseFloat(document.getElementById("bondDiameter").value);
	// Dangerous, assumes that there is a bond for the first atom
	let bondGeometry = particleGroup.children[0].children[1].geometry;
	// This does not work for the box, only for the spheres!
	let scale = diameter / bondGeometry.boundingSphere.radius;
	// console.log(bondGeometry)
	bondGeometry.scale(scale, scale, 1);
};

document.getElementById("resolution").onchange = function () {
	let resolution = parseInt(document.getElementById("resolution").value);
	config["resolution"] = resolution;
	build_scene(animation_frame);
};

document.getElementById("resolution").oninput = function () {
	let resolution = parseInt(document.getElementById("resolution").value);
	document.getElementById("resolutionLabel").innerHTML = "Resolution: " + resolution;
};

o_materialSelect.onchange = function () {
	build_scene(animation_frame);
}

o_wireframe.onchange = function () {
	build_scene(animation_frame);
}

o_spotLightIntensity.oninput = function () {
	spotLight.intensity = o_spotLightIntensity.value;
	document.getElementById("spotLightIntensityLabel").innerHTML = "Spot light intensity: " + o_spotLightIntensity.value;
}

o_hemisphereLightIntensity.oninput = function () {
	hemisphereLight.intensity = o_hemisphereLightIntensity.value;
	document.getElementById("hemisphereLightIntensityLabel").innerHTML = "Hemisphere light intensity: " + o_hemisphereLightIntensity.value;

}

// o_help_btn.onmouseover = function () {
// 	div_help_container.style.display = "block";
// }

// o_help_btn.onmouseout = function () {
// 	div_help_container.style.display = "none";
// }

// o_add_btn.onclick = function () {
// 	document.getElementById("add_class").style.display = "block";
// }


/**
 * Helper function, move later
 * @param {} name 
 * @param {*} checked 
 * @returns 
 */
function createRadioElement(id, properties) {

	var radioFragment = document.createElement('div');
	radioFragment.classList.add("mb-3");
	radioFragment.classList.add("collapse", "show", "scene-modifier");
	radioFragment.id = "sceneModifier_" + id;

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
			fetch("set_update_function_parameter", {
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

addSceneModifier.onchange = function () {
	console.log(this.value);
	if (this.value == "add") {
		addModifierModal.show();
	} else {
		fetch("/select_update_function/" + this.value)
	}

	let domElements = document.getElementsByClassName("scene-modifier");

	[...domElements].forEach(element => {
		let bs_collapse = new bootstrap.Collapse(element, {
			toggle: false
		});
		if (element.id == "sceneModifier_" + this.value) {
			bs_collapse.show();
		} else {
			bs_collapse.hide();
		}
	});
};

document.getElementById("addSceneModifierImportBtn").onclick = function () {
	let function_id = document.getElementById("addSceneModifierImport").value;
	fetch("add_update_function", {
		"method": "POST",
		"headers": { "Content-Type": "application/json" },
		"body": JSON.stringify(function_id),
	}).then(response => response.json()).then(function (response_json) {
		// if not null alert
		if ("error" in response_json) {
			// TODO check if method is already loaded
			alert(response_json["error"]);
			stepError(response_json["error"]);
		} else {
			if (document.getElementById("sceneModifier_" + response_json["title"]) != null) {
				alert("Function already loaded");
				stepError("Function already loaded");
			};

			addModifierModal.hide();

			console.log(response_json);
			load_config();
		}
		return response_json;
	}).then(function (response_json) {
		console.log(response_json);
		let modifier = document.createElement("option");
		modifier.value = response_json["title"];
		modifier.innerHTML = response_json["title"];
		addSceneModifier.appendChild(modifier);
		addSceneModifier.value = response_json["title"];
		return response_json;
	}).then(function (response_json) {
		let sceneModifierSettings = document.getElementById("sceneModifierSettings");
		// sceneModifierSettings.innerHTML = "";
		// TODO make them invisible and only the select one displayed / collapse / none ?
		sceneModifierSettings.appendChild(createRadioElement(response_json["title"], response_json["properties"]));
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
	if (event.isComposing || event.key === "i") {
		printIndices(camera);
	};
});

window.addEventListener("keyup", (event) => {
	if (event.isComposing || event.key === "i") {
		let ids = document.getElementsByClassName("particle-id")
		while (ids.length > 0) {
			ids[0].parentNode.removeChild(ids[0]);
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
		}).then((response) => {
			frames.length = animation_frame + 1;

			getAnimationFrames();
		});

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

	if (frames[animation_frame].length != particleGroup.children.length) {
		// we need to update the scene
		build_scene(animation_frame);
		scene_building = true;
		return; // we need to wait for the scene to be updated
	}

	div_info.innerHTML = "Frame (" + animation_frame + "/" + (frames.length - 1) + ")";

	if (animation_frame != displayed_frame) {
		displayed_frame = animation_frame;
		updateParticlePositions(frames[animation_frame]);
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
		controls.target = getAtomById(selected_ids[0]).position.clone();
	}
}

/**
 * Dynamic indices
 * 
 */


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
