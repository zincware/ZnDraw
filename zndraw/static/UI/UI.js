import { materials } from "../World/components/particles.js";
import { createElementFromSchema } from "./schemaforms.js";

const addModifierModal = new bootstrap.Modal(
  document.getElementById("addModifierModal"),
);

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

async function addSceneModifierOption(function_id) {
  await fetch("add_update_function", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(function_id),
  })
    .then((response) => response.json())
    .then(function (response_json) {
      // if not null alert
      if ("error" in response_json) {
        // TODO check if method is already loaded
        alert(response_json["error"]);
        stepError(response_json["error"]);
      } else {
        if (
          document.getElementById("scene-modifier_" + response_json["title"]) !=
          null
        ) {
          alert("Function already loaded");
          stepError("Function already loaded");
        }
        addModifierModal.hide();
      }
      return response_json;
    })
    .then(function (response_json) {
      let modifier = document.createElement("option");
      modifier.value = response_json["title"];
      modifier.innerHTML = response_json["title"];
      document.getElementById("addSceneModifier").appendChild(modifier);
      return response_json;
    })
    .then(function (response_json) {
      let sceneModifierSettings = document.getElementById(
        "sceneModifierSettings",
      );
      sceneModifierSettings.appendChild(
        createElementFromSchema(response_json, "scene-modifier"),
      );
      document.getElementById("addSceneModifier").value =
        response_json["title"];
    });
}

// load analysis methods from config
async function loadSceneModifier(config) {
  // iterate DATA.config.analysis_methods and add them to the select
  for (let modify_function of config.config.modify_functions) {
    try {
      await addSceneModifierOption(modify_function);
    } catch (error) {
      console.log(error);
    }
  }
  document.getElementById("addSceneModifier").value = "";
  document
    .getElementById("addSceneModifier")
    .dispatchEvent(new Event("change"));

  document.getElementById("addSceneModifier").onchange = function () {
    console.log(this.value);
    if (this.value == "add") {
      addModifierModal.show();
    }

    let domElements = document.getElementsByClassName("scene-modifier");

    [...domElements].forEach((element) => {
      let bs_collapse = new bootstrap.Collapse(element, {
        toggle: false,
      });
      if (element.id == "scene-modifier_" + this.value) {
        bs_collapse.show();
      } else {
        bs_collapse.hide();
      }
    });
  };

  document.getElementById("sceneModifierBtn").onclick = function () {
    // div_info.innerHTML = "Processing...";

    let form = document.getElementById(
      "scene-modifier_" + document.getElementById("addSceneModifier").value,
    );
    let modifier_kwargs = {};
    Array.from(form.elements).forEach((input) => {
      modifier_kwargs[input.dataset.key] = input.value;
    });

    fetch("update", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        selected_ids: config.selected,
        step: config.step,
        modifier: document.getElementById("addSceneModifier").value,
        modifier_kwargs: modifier_kwargs,
      }),
    }).then(() => {
      config.deleteCache();
    });
  };
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
  loadSceneModifier(config);

  // disable loading spinner by making it invisible
  const loadingElem = document.getElementById("atom-spinner");
  loadingElem.style.display = "none";
}
