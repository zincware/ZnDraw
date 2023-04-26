// Convert schema to bootstrap form


export function createElementFromSchema(schema, clsName) {

	let modifierCanvas = document.createElement('form');
	// modifierCanvas.classList.add("mb-3");
	modifierCanvas.classList.add("collapse", "show", clsName, "border", "border-primary", "rounded", "p-3");
	modifierCanvas.id = clsName+ "_" + schema.title;

	console.log("Adding modifier: " + schema.title);

	Object.values(schema.properties).forEach((item) => {
		let controller = document.createElement('input');
		controller.dataset.key = item["title"];
		if (item["type"] == "integer") {
			controller.type = "range";
			controller.classList.add("form-range");
			controller.step = 1;
		} else if (item["type"] == "number") {
			controller.type = "range";
			controller.classList.add("form-range");
			controller.step = 0.1;
		} else if (item["type"] == "text") {
			controller.type = "text";
			controller.classList.add("form-control");
		} else {
			console.log("Unknown type: " + item["type"]);
		}
		controller.value = item["default"];
		controller.id = schema.title + "_" + item["title"];

		if ("minimum" in item) {
			controller.min = item["minimum"];
		}
		if ("maximum" in item) {
			controller.max = item["maximum"];
		}

		let controller_label = document.createElement('label');
		controller_label.classList.add("form-label");
		controller_label.setAttribute("for", controller.id);
		controller_label.innerHTML = item["title"] + ": " + controller.value;

		let function_container = document.createElement('div');
		// function_container.classList.add("mb-1");
		function_container.appendChild(controller_label);
		function_container.appendChild(controller);

		modifierCanvas.appendChild(function_container);

		controller.oninput = function () {
			controller_label.innerHTML = item["title"] + ": " + controller.value;
		}
	});

	return modifierCanvas;
}
