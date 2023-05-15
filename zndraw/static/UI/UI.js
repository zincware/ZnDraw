import { materials } from "../World/components/particles.js";

function update_materials(config) {
  const o_materialSelect = document.getElementById("materialSelect");

  for (const material in materials) {
    const option = document.createElement("option");
    option.text = material;
    option.value = material;
    o_materialSelect.appendChild(option);
  }
  o_materialSelect.value = config.config.material;

  o_materialSelect.onchange = function () {
    config.update({ material: o_materialSelect.value });
  };

  document.getElementById("wireframe").onchange = function () {
    config.update({ material_wireframe: this.checked });
  };
}

function updateFPS(config) {
  document.getElementById("max_fps").onchange = function () {
    config.update({ max_fps: this.value });
  };
}

export function setUIEvents(config) {
  update_materials(config);
  updateFPS(config);

  // disable loading spinner by making it invisible
  const loadingElem = document.getElementById("atom-spinner");
  loadingElem.style.display = "none";
}
