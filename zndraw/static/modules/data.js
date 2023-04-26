// Handle Data Loading

export const frames = { "position": [], "force": [], length: 0, "box": [] };
let data_loading = false;
export let rebuild_cache = {};


export function spliceFrames(stop) {
	let start = 0;
	frames.position.splice(start, stop);
	frames.force.splice(start, stop);
	frames.box.splice(start, stop);
	frames.length = frames.position.length;
}

export let config = {};

let _config_update = false;

const frames_per_post = 100;

// Move config to this module?
export async function load_config() {
	if (_config_update) {
		return;
	}
	_config_update = true;
	config = await (await fetch("config")).json();
	console.log(config)
	_config_update = false;
}

export async function getAnimationFrames() {
	if (data_loading) {
		console.log("Animation read already in progress");
		return;
	}
	data_loading = true;
	console.log("Loading animation frames");

	if (frames.position.length < 2) {
		fetch("load");
	}

	let step = frames.position.length;
	while (true) { // TODO make non-stop loading optional

		let obj = await fetch("data", {
			"method": "POST",
			"headers": { "Content-Type": "application/json" },
			"body": JSON.stringify({ "start": step, "stop": frames_per_post + step }), // TODO keep configurable
		}).then(response => response.json()).then(function (response_json) {
			load_config();
			return response_json;
		});

		console.log("Read " + step + "-" + (frames_per_post + step) + " frames");
		if (obj["position"].length === 0) {
			if (!config["continuous_loading"]){
				console.log("Animation read finished");
				break;
			} else {
				console.log("Continuous Loading: waiting for new frames");
			}
		}
		frames.position = frames.position.concat(obj["position"]);
		frames.box = frames.box.concat(obj["box"]);
		// frames.force = frames.force.concat(obj["force"]);
		frames.length = frames.position.length;

		step = frames.length;
	}
	data_loading = false;
}

export async function getRebuildCache(step) {
	let arrayOfResponses = [];
	if (rebuild_cache.hasOwnProperty(step)) {
		console.log("Using cached scene");
		arrayOfResponses = rebuild_cache[step];
	} else {
		// this is faster then doing it one by one

		arrayOfResponses = await (await fetch("graph", {
			"method": "POST",
			"headers": { "Content-Type": "application/json" },
			"body": JSON.stringify(step)
		})).json();
		rebuild_cache[step] = arrayOfResponses;
	}
	await load_config();
	if( (step === 0) && ( frames.position.length === 0 )){
		let positions = [];
		for (let i = 0; i < arrayOfResponses["nodes"].length; i++) {
			positions.push(arrayOfResponses["nodes"][i]["position"]);
		}
		frames.position = [positions];
		frames.box = [arrayOfResponses["box"]];
		frames.length = frames.position.length;
	}
	return arrayOfResponses;
}

export function resetAnimationFrames() {
	console.log("Resetting animation frames")
	frames.position = [];
	frames.force = [];
	frames.length = 0;
	frames.box = [];
	rebuild_cache = {};
	getRebuildCache(0);
	getAnimationFrames();
}
