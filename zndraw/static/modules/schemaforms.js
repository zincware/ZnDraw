// Convert schema to bootstrap form


export function createElementFromSchema(schema, clsName) {

	let modifierCanvas = document.createElement('form');
	// modifierCanvas.classList.add("mb-3");
	modifierCanvas.classList.add("collapse", "show", clsName, "border", "border-primary", "rounded", "p-3");
	modifierCanvas.id = clsName+ "_" + schema.title;

	console.log("Adding modifier: " + schema.title);
	console.log(schema);

	for (const [key, item] of Object.entries(schema.properties)) {
		let controller = document.createElement('input');
		controller.dataset.key = key;
		if (item["type"] == "integer") {
			controller.type = "range";
			controller.classList.add("form-range");
			controller.step = 1;
		} else if (item["type"] == "number") {
			controller.type = "range";
			controller.classList.add("form-range");
			controller.step = 0.1;
		} else if (item["type"] == "boolean") {
			controller.type = "checkbox";
			controller.classList.add("form-check-input");
			controller.checked = item["default"];
			controller.onchange = function () {
				this.value = this.checked;
			}
			// controller.disabled = true;
		} else if (["test", "string"].includes(item["type"])) {
			// check if enum is available to create a select
			if ("enum" in item) {
				// we overwrite the controller with a select!
				controller = document.createElement('select');
				controller.classList.add("form-select");
				controller.dataset.key = key;
				item["enum"].forEach((enum_item) => {
					let option = document.createElement('option');
					option.value = enum_item;
					option.innerHTML = enum_item;
					controller.appendChild(option);
				});
			} else {
				controller.type = "text";
				controller.classList.add("form-control");
			}
		} else {
			console.log("Unknown type: " + item["type"]);
		}
		controller.value = item["default"];
		// controller.id = schema.title + "_" + item["title"];

		if ("minimum" in item) {
			controller.min = item["minimum"];
		}
		if ("maximum" in item) {
			controller.max = item["maximum"];
		}

		let controller_label = document.createElement('label');
		controller_label.classList.add("form-label");
		controller_label.setAttribute("for", controller.id);
		if (["text", "boolean"].includes(item["type"])) {
			controller_label.innerHTML = item["title"];
		} else {
			controller_label.innerHTML = item["title"] + ": " + controller.value;
			controller.oninput = function () {
				controller_label.innerHTML = item["title"] + ": " + controller.value;
			}
		}

		let function_container = document.createElement('div');
		// function_container.classList.add("mb-1");
		function_container.appendChild(controller_label);
		function_container.appendChild(controller);

		modifierCanvas.appendChild(function_container);
	};

	return modifierCanvas;
}
