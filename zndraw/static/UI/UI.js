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

function update_resolution(config) {
  const o_resolution = document.getElementById("resolution");
  o_resolution.value = config.config.resolution;

  o_resolution.onchange = function () {
    config.update({ resolution: parseInt(o_resolution.value) }).then(() => {
      config.rebuild();
    });
    document.getElementById("resolutionLabel").innerHTML =
      "Resolution: " + this.value;
  };
}

function update_sphere_radius(config) {
  const o_sphere_radius = document.getElementById("sphereRadius");
  o_sphere_radius.value = config.config.sphere_radius;

  o_sphere_radius.onchange = function () {
    config.update({ sphere_size: parseFloat(o_sphere_radius.value) });
    document.getElementById("sphereRadiusLabel").innerHTML =
      "Sphere radius: " + this.value;
  };
}

function update_bond_radius(config) {
  const o_bond_radius = document.getElementById("bondDiameter");
  o_bond_radius.value = config.config.bond_radius;

  o_bond_radius.onchange = function () {
    config.update({ bond_size: parseFloat(o_bond_radius.value) }).then(() => {
      config.rebuild();
    });
    document.getElementById("bondDiameterLabel").innerHTML =
      "Bond diameter: " + this.value;
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

    if (event.isComposing || event.key === "ArrowRight") {
      config.play = false;
      config.set_step(
        Math.min(config.config.total_frames - 1, config.step + 1),
      );
    }
    if (event.isComposing || event.key === "ArrowLeft") {
      config.play = false;
      config.set_step(Math.max(0, config.step - 1));
    }
    if (event.isComposing || event.key === "ArrowUp") {
      config.play = false;
      config.set_step(
        parseInt(
          Math.min(
            config.config.total_frames - 1,
            config.step + config.config.total_frames / 10,
          ),
        ),
      );
    }
    if (event.isComposing || event.key === "ArrowDown") {
      config.play = false;
      config.set_step(
        parseInt(Math.max(0, config.step - config.config.total_frames / 10)),
      );
    }
  });
}

function attachKeyPressed(config) {
  window.addEventListener("keydown", (event) => {
    config.pressed_keys[event.key] = true;
  });
  window.addEventListener("keyup", (event) => {
    config.pressed_keys[event.key] = false;
  });
}

export function setUIEvents(config) {
  update_materials(config);
  updateFPS(config);
  setupSlider(config);
  setupPlayPause(config);
  attachKeyPressed(config);
  update_resolution(config);
  update_sphere_radius(config);
  update_bond_radius(config);

  // disable loading spinner by making it invisible
  const loadingElem = document.getElementById("atom-spinner");
  loadingElem.style.display = "none";
}
