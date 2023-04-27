import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import * as PARTICLES from './modules/particles.js';
import * as DRAW from './modules/draw.js';
import * as DATA from './modules/data.js';
import { keydown, keyconfig } from './modules/keypress.js';
import { createElementFromSchema } from './modules/schemaforms.js';
// THREE.Cache.enabled = true;


// some global variables
let selected_ids = [];
let animation_frame = 0;
let displayed_frame = 0;
let scene_building = false;
let animation_running = false;
let fps = [];


/**
 * DOM variables
 */

const div_info = document.getElementById('info');
const div_loading = document.getElementById('atom-spinner');
const div_greyOut = document.getElementById('greyOut');
const div_lst_selected_ids = document.getElementById('lst_selected_ids');
const div_FPS = document.getElementById('FPS');
const div_n_particles = document.getElementById('n_particles');
const div_n_bonds = document.getElementById('n_bonds');

const o_autoRestart = document.getElementById('autoRestart');
const o_animate = document.getElementById('animate');
const o_reset_selection = document.getElementById('reset_selection');
const o_reset = document.getElementById('reset');
const o_max_fps = document.getElementById('max_fps');
const o_materialSelect = document.getElementById('materialSelect');
const o_wireframe = document.getElementById('wireframe');
const o_spotLightIntensity = document.getElementById('spotLightIntensity');
const o_hemisphereLightIntensity = document.getElementById('hemisphereLightIntensity');

const addModifierModal = new bootstrap.Modal(document.getElementById("addModifierModal"));
const addAnalysisModal = new bootstrap.Modal(document.getElementById("addAnalysisModal"));
const addSceneModifier = document.getElementById("addSceneModifier");


// Helper Functions

await DATA.load_config();
// THREE JS Setup
const scene = new THREE.Scene();

function getCamera() {
	let width = window.innerWidth;
	let height = window.innerHeight;

	if (DATA.config.camera == "PerspectiveCamera") {
		return new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
	} else {
		return new THREE.OrthographicCamera(width / - 2, width / 2, height / 2, height / - 2, 1, 1000);
	}
}

const camera = getCamera();
camera.position.z = 50;

const renderer = new THREE.WebGLRenderer({ antialias: DATA.config["antialias"], alpha: true });

const hemisphereLight = new THREE.HemisphereLight(0xffffff, 0x777777, 0.1);
const spotLight = new THREE.SpotLight(0xffffff, 0, 0, Math.PI / 2);
const cameraLight = new THREE.PointLight(0xffffff, 1.0);

camera.add(cameraLight);

renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);


spotLight.position.set(0, 0, 100);
scene.add(spotLight);
scene.add(camera); // add the camera to the scene because of the cameraLight

scene.add(hemisphereLight);

async function update_materials() {
	for (const material in PARTICLES.materials) {
		const option = document.createElement("option");
		option.text = material;
		option.value = material;
		o_materialSelect.appendChild(option);
	}
	o_materialSelect.value = DATA.config["material"];
}

update_materials();



// Setup Scene


async function build_scene(step) {
	if (scene_building) {
		return;
	}
	scene_building = true;
	console.log("Updating scene");

	div_loading.style.visibility = 'visible';
	div_greyOut.style.visibility = 'visible';

	animation_frame = step;

	await DATA.getRebuildCache(step).then(function (arrayOfResponses) {
		PARTICLES.drawAtoms(arrayOfResponses["nodes"], arrayOfResponses["edges"], DATA.config, scene);
		selected_ids = [];
	}).then(update_selection).then(function () {
		scene_building = false;
		div_n_particles.innerHTML = PARTICLES.particleGroup.children.length;
		div_n_bonds.innerHTML = PARTICLES.countBonds();
		div_greyOut.style.visibility = 'hidden';
		div_loading.style.visibility = 'hidden';
	});
}

build_scene(0);

// interactions

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
	const intersects = raycaster.intersectObjects(PARTICLES.particleGroup.children, true);

	if (intersects.length == 0) {
		return;
	}

	// for (let i = 0; i < intersects.length; i++) {
	let mesh = intersects[0].object;
	if (!keydown["shift"] && !keyconfig.multiselect) {
		if (selected_ids.includes(mesh.userData["id"])) {
			selected_ids = [];
		} else {
			selected_ids = [mesh.userData["id"]];
		}
	} else {
		if (!selected_ids.includes(mesh.userData["id"])) {
			selected_ids.push(mesh.userData["id"]);
		} else {
			selected_ids.splice(selected_ids.indexOf(mesh.userData["id"]), 1);
		}
	}
	// update colors here for better performance
	// TODO chain them
	await update_color_of_ids(selected_ids);
	await update_selection();
}

/**
 * We update the color of every atom in the scene
 * @param {list[int]} ids, the selected ids
 * @returns 
 */

async function update_color_of_ids(ids) {
	PARTICLES.particleGroup.children.forEach(function (atom_grp) {
		let mesh = atom_grp.children[0];
		if (ids.includes(mesh.userData["id"])) {
			let material = PARTICLES.speciesMaterial(document.getElementById('materialSelect').value, 0xffa500, document.getElementById('wireframe').checked);
			atom_grp.children.forEach(function (item) {
				item.material = material;
			});
		} else {
			let material = PARTICLES.speciesMaterial(document.getElementById('materialSelect').value, mesh.userData["color"], document.getElementById('wireframe').checked);
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


if (DATA.config["restart_animation"] === true) {
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
	DATA.getAnimationFrames();
}


o_reset_selection.onclick = function () {
	selected_ids = [];
	update_color_of_ids(selected_ids);
	update_selection();
}


o_reset.onclick = function () {
	DATA.load_config();
	animation_frame = 0;
	build_scene(0);
	selected_ids = [];
	update_selection();
}

document.getElementById("sphereRadius").onchange = function () {
	let radius = parseFloat(document.getElementById("sphereRadius").value);
	let particleGeometry = PARTICLES.particleGroup.children[0].children[0].geometry;
	let scale = radius / particleGeometry.boundingSphere.radius;
	particleGeometry.scale(scale, scale, scale);
	fetch("config", {
		"method": "POST",
		"headers": { "Content-Type": "application/json" },
		"body": JSON.stringify({ "sphere_size": radius }),
	});

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
	fetch("config", {
		"method": "POST",
		"headers": { "Content-Type": "application/json" },
		"body": JSON.stringify({ "bond_size": diameter }),
	});
	// currently we can't scale correctly so we reset everything
	build_scene(animation_frame);
};

document.getElementById("resolution").onchange = function () {
	let resolution = parseInt(document.getElementById("resolution").value);
	fetch("config", {
		"method": "POST",
		"headers": { "Content-Type": "application/json" },
		"body": JSON.stringify({ "resolution": resolution }),
	});
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

document.getElementById("cameraLightIntensity").oninput = function () {
	let intensity = parseFloat(document.getElementById("cameraLightIntensity").value);
	cameraLight.intensity = intensity;
	document.getElementById("cameraLightIntensityLabel").innerHTML = "Camera light intensity: " + intensity;
}



addSceneModifier.onchange = function () {
	console.log(this.value);
	if (this.value == "add") {
		addModifierModal.show();
	}

	let domElements = document.getElementsByClassName("scene-modifier");

	[...domElements].forEach(element => {
		let bs_collapse = new bootstrap.Collapse(element, {
			toggle: false
		});
		if (element.id == "scene-modifier_" + this.value) {
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
			if (document.getElementById("scene-modifier_" + response_json["title"]) != null) {
				alert("Function already loaded");
				stepError("Function already loaded");
			};
			addModifierModal.hide();
			DATA.load_config();
		}
		return response_json;
	}).then(function (response_json) {
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
		sceneModifierSettings.appendChild(createElementFromSchema(response_json, "scene-modifier"));
	});
}


document.getElementById("addAnalysis").onchange = function () {
	console.log(this.value);
	if (this.value == "add") {
		addAnalysisModal.show();
	}

	let domElements = document.getElementsByClassName("scene-analysis");

	[...domElements].forEach(element => {
		let bs_collapse = new bootstrap.Collapse(element, {
			toggle: false
		});
		if (element.id == "scene-analysis_" + this.value) {
			bs_collapse.show();
		} else {
			bs_collapse.hide();
		}
	});
};

document.getElementById("addAnalysisImportBtn").onclick = function () {
	let function_id = document.getElementById("addAnalysisImport").value;
	fetch("add_analysis", {
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
			if (document.getElementById("scene-analysis_" + response_json["title"]) != null) {
				alert("Function already loaded");
				stepError("Function already loaded");
			};
			addAnalysisModal.hide();
			DATA.load_config();
		}
		return response_json;
	}).then(function (response_json) {
		let modifier = document.createElement("option");
		modifier.value = response_json["title"];
		modifier.innerHTML = response_json["title"];
		document.getElementById("addAnalysis").appendChild(modifier);
		document.getElementById("addAnalysis").value = response_json["title"];
		return response_json;
	}).then(function (response_json) {
		let analysisSettings = document.getElementById("analysisSettings");
		analysisSettings.appendChild(createElementFromSchema(response_json, "scene-analysis"));
	});
}

window.addEventListener("keydown", (event) => {
	if (event.isComposing || event.key === " ") {
		if (document.activeElement == document.body) {
			if (animation_running && animation_frame == DATA.frames.length - 1) {
				animation_frame = 0;
			} else {
				animation_running = !animation_running;
			}
		}
	}
	if (event.isComposing || event.key === "ArrowLeft") {
		if (document.activeElement == document.body) {
			animation_running = false;
			animation_frame = Math.max(0, animation_frame - 1);
		}
	}
	if (event.isComposing || event.key === "ArrowRight") {
		if (document.activeElement == document.body) {
		animation_running = false;
		animation_frame = Math.min(DATA.frames.length - 1, animation_frame + 1);
		}
	}
	if (event.isComposing || event.key === "ArrowUp") {
		if (document.activeElement == document.body) {
			animation_running = false;
			animation_frame = parseInt(Math.min(DATA.frames.length - 1, animation_frame + (DATA.frames.length / 10)));
		}
	}
	if (event.isComposing || event.key === "ArrowDown") {
		if (document.activeElement == document.body) {
			animation_running = false;
			animation_frame = parseInt(Math.max(0, animation_frame - (DATA.frames.length / 10)));
		}
	}
	if (event.isComposing || event.key === "q") {
		DATA.resetAnimationFrames();
		animation_frame = 0;
		document.getElementById("frame-slider").value = 0;
	}
	if (event.isComposing || event.key === "i") {
		if (document.activeElement == document.body) {
			PARTICLES.printIndices(camera);
		}
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

document.getElementById("sceneModifierBtn").onclick = function () {
	div_info.innerHTML = "Processing...";

	let form = document.getElementById("scene-modifier_" + addSceneModifier.value);
	let modifier_kwargs = {}
	Array.from(form.elements).forEach((input) => {
		modifier_kwargs[input.dataset.key] = input.value;
	});

	fetch("update", {
		"method": "POST",
		"headers": { "Content-Type": "application/json" },
		"body": JSON.stringify({ "selected_ids": selected_ids, "step": animation_frame, "points": DRAW.positions, "modifier": addSceneModifier.value, "modifier_kwargs": modifier_kwargs }),
	}).then(function (response) {
		DATA.resetAnimationFrames(); // use DATA.spliceFrames(animation_frame + 1); ?
	}).then(DATA.getAnimationFrames);
}



// Drawing

DRAW.init(camera, renderer, scene, controls)

document.getElementById("drawAddAnchor").onclick = function () {
	let { position } = PARTICLES.particleGroup.children[selected_ids[0]].children[0];
	DRAW.createAnchorPoint(position);
};

document.getElementById("drawRemoveLine").onclick = function () {
	DRAW.reset();
};

document.getElementById("analyseBtn").onclick = function () {
	let form = document.getElementById("scene-analysis_" + document.getElementById("addAnalysis").value);
	let modifier_kwargs = {}
	Array.from(form.elements).forEach((input) => {
		modifier_kwargs[input.dataset.key.toLowerCase()] = input.value;
	});

	fetch("analyse", {
		"method": "POST",
		"headers": { "Content-Type": "application/json" },
		"body": JSON.stringify({ "selected_ids": selected_ids, "step": animation_frame, "points": DRAW.positions, "modifier": document.getElementById("addAnalysis").value, "modifier_kwargs": modifier_kwargs }),
	}).then((response) => response.json()).then(function (response_json) {
		Plotly.newPlot('analysePlot', response_json);
		document.getElementById("analysePlot").on('plotly_click', function (data) {
			console.log(data);
			animation_frame = data.points[0].pointIndex;
		});
	});
};

document.getElementById("showBox").onclick = function () {
	if (this.checked) {
		scene.add(PARTICLES.createBox(DATA.frames.box[animation_frame]));
	} else {
		scene.remove(PARTICLES.box);
	}
};

document.getElementById("continuousLoading").onclick = function () {
	if (this.checked) {
		fetch("config",
			{
				"method": "POST",
				"headers": { "Content-Type": "application/json" },
				"body": JSON.stringify({ "continuous_loading": true }),
			}).then(DATA.load_config).then(DATA.getAnimationFrames);
	} else {
		fetch("config",
			{
				"method": "POST",
				"headers": { "Content-Type": "application/json" },
				"body": JSON.stringify({ "continuous_loading": false }),
			}).then(DATA.load_config);
	}
};

document.getElementById("frame-slider").oninput = function () {
	if (this.value > DATA.frames.length - 1) {
		this.value = DATA.frames.length - 1;
	}
	animation_frame = parseInt(this.value);
};

let move_atoms_clock = new THREE.Clock();


function move_atoms() {
	if (DATA.frames.length < animation_frame) {
		return;
	}
	if (move_atoms_clock.getElapsedTime() < (1 / o_max_fps.value)) {
		return;
	}
	if (scene_building === true) {
		return;
	}

	document.getElementById("frame-slider").value = animation_frame;
	document.getElementById("frame-slider").max = DATA.config["total_frames"];

	try {
		if (DATA.frames.position[animation_frame].length != PARTICLES.particleGroup.children.length) {
			// we need to update the scene
			console.log("Particle count changed from " + DATA.frames.position[animation_frame].length + " to " + PARTICLES.particleGroup.children.length + ", rebuilding scene");
			build_scene(animation_frame);
			return; // we need to wait for the scene to be updated
		}
	} catch (e) {
		console.log("Scene not ready yet");
		return; // scene is not ready yet
	}

	div_info.innerHTML = "Frame (" + animation_frame + "/" + (DATA.frames.length - 1) + ")";

	if (animation_frame != displayed_frame) {
		displayed_frame = animation_frame;
		PARTICLES.updateParticlePositions(DATA.frames.position[animation_frame]);
		// if (PARTICLES.boxGeometry !== undefined) {
		// 	PARTICLES.updateBox(DATA.frames.box[animation_frame]);
		// }
		// PARTICLES.updateArrows(DATA.frames.force[animation_frame]);
	}

	fps.push(move_atoms_clock.getElapsedTime());

	if (fps.length > 10) {
		fps.shift();
	}
	// if (!data_loading) {
	// 	DATA.getAnimationFrames();
	// }

	div_FPS.innerHTML = (1 / (fps.reduce((a, b) => a + b, 0) / fps.length)).toFixed(2);
	move_atoms_clock.start();
	if (animation_running === true) {
		if (animation_frame < DATA.frames.length - 1) {
			animation_frame += 1;
		} else if (o_autoRestart.checked === true) {
			animation_frame = 0;
		}
	}
}

function centerCamera() {
	if (selected_ids.length === 0) {
		controls.target = PARTICLES.getAtomsCenter([...Array(PARTICLES.particleGroup.children.length).keys()]);
	} else {
		controls.target = PARTICLES.getAtomsCenter(selected_ids);
	}
}

div_info.innerHTML = "Reading file...";
DATA.getAnimationFrames();

// scene.add(PARTICLES.arrowGroup);
function animate() {

	renderer.render(scene, camera);
	controls.update();

	move_atoms();
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
