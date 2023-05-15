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

function setupSlider(config) {
  const slider = document.getElementById("frame-slider");
  config.set_step_callbacks.push((config) => {
    slider.value = config.step;
    slider.max = config.config.total_frames;
  });

  slider.oninput = function () {
    config.set_step(parseInt(this.value));
  };
}

function setupPlayPause(config) {
  window.addEventListener("keydown", (event) => {
    if (event.isComposing || event.key === " ") {
      config.play = !config.play;
      if (config.play) {
        console.log("play");
      } else {
        console.log("pause");
      }
    }
  });
}

export function setUIEvents(config) {
  update_materials(config);
  updateFPS(config);
  setupSlider(config);
  setupPlayPause(config);

  // disable loading spinner by making it invisible
  const loadingElem = document.getElementById("atom-spinner");
  loadingElem.style.display = "none";
}
