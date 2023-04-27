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

document.getElementById("multiselect").onchange = function () {
	keyconfig.multiselect = this.checked;
}



// Rescale offcanvas by dragging
let active_offcanvas_border;
const offcanvas_borders = document.getElementsByClassName("offcanvas-border");


function resize_offcanvas(e) {
	if (e.clientX < 200) {
		return;
	}
	active_offcanvas_border.parentNode.style.width = e.clientX + "px";
	active_offcanvas_border.style.left = e.clientX + "px";
}

for (let i = 0; i < offcanvas_borders.length; i++) {
	offcanvas_borders[i].onpointerdown = function (e) {
		console.log(this);
		active_offcanvas_border = this;
		document.addEventListener("pointermove", resize_offcanvas);
	}
}

document.addEventListener("pointerup", function (e) {
	document.removeEventListener("pointermove", resize_offcanvas);
});

