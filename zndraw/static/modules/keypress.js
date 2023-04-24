export const keydown = { "shift": false, "ctrl": false, "alt": false, "c": false, "l": false };
export const keyconfig = { "multiselect": false }



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

document.getElementById("multiselect").onchange = function() {
	keyconfig.multiselect = this.checked;
}