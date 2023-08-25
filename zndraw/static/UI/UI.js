import { materials } from '../World/components/particles.js';

const addModifierModal = new bootstrap.Modal(
  document.getElementById('addModifierModal'),
);

const addAnalysisModal = new bootstrap.Modal(
  document.getElementById('addAnalysisModal'),
);

function update_materials(config) {
  const o_materialSelect = document.getElementById('materialSelect');

  for (const material in materials) {
    const option = document.createElement('option');
    option.text = material;
    option.value = material;
    o_materialSelect.appendChild(option);
  }
  o_materialSelect.value = config.config.material;

  o_materialSelect.onchange = function () {
    config.update({ material: o_materialSelect.value });
  };

  document.getElementById('wireframe').onchange = function () {
    config.update({ material_wireframe: this.checked });
  };
}

function update_resolution(config, world) {
  const o_resolution = document.getElementById('resolution');
  o_resolution.value = config.config.resolution;

  o_resolution.onchange = function () {
    config.update({ resolution: parseInt(o_resolution.value) }).then(() => {
      world.rebuild();
    });
    document.getElementById(
      'resolutionLabel',
    ).innerHTML = `Resolution: ${this.value}`;
  };
}

function update_sphere_radius(config) {
  const o_sphere_radius = document.getElementById('sphereRadius');
  o_sphere_radius.value = config.config.sphere_radius;

  o_sphere_radius.onchange = function () {
    config.update({ sphere_size: parseFloat(o_sphere_radius.value) });
    document.getElementById(
      'sphereRadiusLabel',
    ).innerHTML = `Sphere radius: ${this.value}`;
  };
}

function update_bond_radius(config, world) {
  const o_bond_radius = document.getElementById('bondDiameter');
  o_bond_radius.value = config.config.bond_radius;

  o_bond_radius.onchange = function () {
    config.update({ bond_size: parseFloat(o_bond_radius.value) }).then(() => {
      world.rebuild();
    });
    document.getElementById(
      'bondDiameterLabel',
    ).innerHTML = `Bond diameter: ${this.value}`;
  };
}

function updateFPS(config) {
  document.getElementById('max_fps').onchange = function () {
    config.update({ max_fps: this.value });
  };
}

function setupPlayPause(world) {
  window.addEventListener('keydown', (event) => {
    if (event.isComposing || event.key === 'ArrowRight') {
      world.setStep(world.step + 1);
    }
  });
}

// function setupPlayPause(config) {
//   window.addEventListener("keydown", (event) => {
//     if (event.isComposing || event.key === " ") {
//       config.play = !config.play;
//       if (config.play) {
//         console.log("play");
//       } else {
//         console.log("pause");
//       }
//     }

//     if (event.isComposing || event.key === "ArrowRight") {
//       config.play = false;
//       config.set_step(
//         Math.min(config.config.total_frames - 1, config.step + 1),
//       );
//     }
//     if (event.isComposing || event.key === "ArrowLeft") {
//       config.play = false;
//       config.set_step(Math.max(0, config.step - 1));
//     }
//     if (event.isComposing || event.key === "ArrowUp") {
//       config.play = false;
//       config.set_step(
//         parseInt(
//           Math.min(
//             config.config.total_frames - 1,
//             config.step + config.config.total_frames / 10,
//           ),
//         ),
//       );
//     }
//     if (event.isComposing || event.key === "ArrowDown") {
//       config.play = false;
//       config.set_step(
//         parseInt(Math.max(0, config.step - config.config.total_frames / 10)),
//       );
//     }
//     if (event.isComposing || event.key === "End") {
//       config.play = false;
//       config.set_step(config.config.total_frames - 1);
//     }
//   });
// }

function attachKeyPressed(config) {
  window.addEventListener('keydown', (event) => {
    config.pressed_keys[event.key] = true;
  });
  window.addEventListener('keyup', (event) => {
    config.pressed_keys[event.key] = false;
  });
}

async function addSceneModifierOption(function_id) {
  await fetch('add_update_function', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(function_id),
  })
    .then((response) => response.json())
    .then((response_json) => {
      // if not null alert
      if ('error' in response_json) {
        // TODO check if method is already loaded
        alert(response_json.error);
        stepError(response_json.error);
      } else {
        if (
          document.getElementById(`scene-modifier_${response_json.title}`)
          != null
        ) {
          alert('Function already loaded');
          stepError('Function already loaded');
        }
        addModifierModal.hide();
      }
      return response_json;
    })
    .then((response_json) => {
      const modifier = document.createElement('option');
      modifier.value = response_json.title;
      modifier.innerHTML = response_json.title;
      document.getElementById('addSceneModifier').appendChild(modifier);
      return response_json;
    })
    .then((response_json) => {
      const sceneModifierSettings = document.getElementById(
        'sceneModifierSettings',
      );
      sceneModifierSettings.appendChild(
        createElementFromSchema(response_json, 'scene-modifier'),
      );
      document.getElementById('addSceneModifier').value = response_json.title;
    });
}

async function addAnalysisOption(function_id) {
  await fetch('add_analysis', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(function_id),
  })
    .then((response) => response.json())
    .then((response_json) => {
      // if not null alert
      if ('error' in response_json) {
        // TODO check if method is already loaded
        console.log(
          `Adding analysis failed with error: ${response_json.error}`,
        );
        // alert(response_json["error"]);
        // stepError(response_json["error"]);
      } else {
        if (
          document.getElementById(`scene-analysis_${response_json.title}`)
          != null
        ) {
          alert('Function already loaded');
          stepError('Function already loaded');
        }
        addAnalysisModal.hide();
      }
      return response_json;
    })
    .then((response_json) => {
      const modifier = document.createElement('option');
      modifier.value = response_json.title;
      modifier.innerHTML = response_json.title;
      document.getElementById('addAnalysis').appendChild(modifier);
      return response_json;
    })
    .then((response_json) => {
      console.log(response_json);
      const sceneModifierSettings = document.getElementById('analysisSettings');
      sceneModifierSettings.appendChild(
        createElementFromSchema(response_json, 'scene-analysis'),
      );
      document.getElementById('addAnalysis').value = response_json.title;
    });
}

// load analysis methods from config
async function loadSceneModifier(config, world) {
  // iterate DATA.config.analysis_methods and add them to the select
  for (const modify_function of config.config.modify_functions) {
    try {
      await addSceneModifierOption(modify_function);
    } catch (error) {
      console.log(error);
    }
  }
  document.getElementById('addSceneModifier').value = '';
  document
    .getElementById('addSceneModifier')
    .dispatchEvent(new Event('change'));

  document.getElementById('addSceneModifier').onchange = function () {
    console.log(this.value);
    if (this.value == 'add') {
      addModifierModal.show();
    }

    const domElements = document.getElementsByClassName('scene-modifier');

    [...domElements].forEach((element) => {
      const bs_collapse = new bootstrap.Collapse(element, {
        toggle: false,
      });
      if (element.id == `scene-modifier_${this.value}`) {
        bs_collapse.show();
      } else {
        bs_collapse.hide();
      }
    });
  };

  document.getElementById('sceneModifierBtn').onclick = function () {
    // div_info.innerHTML = "Processing...";

    const form = document.getElementById(
      `scene-modifier_${document.getElementById('addSceneModifier').value}`,
    );
    const modifier_kwargs = {};
    Array.from(form.elements).forEach((input) => {
      modifier_kwargs[input.dataset.key] = input.value;
    });

    fetch('update', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        selected_ids: config.selected,
        step: config.step,
        modifier: document.getElementById('addSceneModifier').value,
        modifier_kwargs,
        points: config.draw_vectors,
      }),
    }).then(() => {
      world.deleteCache();
    });
  };
}

async function loadSceneAnalysis(config, world) {
  // iterate DATA.config.analysis_methods and add them to the select
  for (const method of config.config.analysis_functions) {
    try {
      await addAnalysisOption(method);
    } catch (error) {
      console.log(error);
    }
  }
  document.getElementById('addAnalysis').value = '';
  document.getElementById('addAnalysis').dispatchEvent(new Event('change'));

  document.getElementById('addAnalysis').onchange = function () {
    console.log(this.value);
    if (this.value == 'add') {
      addAnalysisModal.show();
    }

    const domElements = document.getElementsByClassName('scene-analysis');
    console.log(domElements);

    [...domElements].forEach((element) => {
      const bs_collapse = new bootstrap.Collapse(element, {
        toggle: false,
      });
      if (element.id == `scene-analysis_${this.value}`) {
        bs_collapse.show();
      } else {
        bs_collapse.hide();
      }
    });
  };

  document.getElementById('analyseBtn').onclick = function () {
    // div_info.innerHTML = "Processing...";

    const form = document.getElementById(
      `scene-analysis_${document.getElementById('addAnalysis').value}`,
    );
    const modifier_kwargs = {};
    Array.from(form.elements).forEach((input) => {
      modifier_kwargs[input.dataset.key] = input.value;
    });

    fetch('analyse', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        selected_ids: config.selected,
        step: config.step,
        points: [],
        modifier: document.getElementById('addAnalysis').value,
        modifier_kwargs,
      }),
    })
      .then((response) => response.json())
      .then((response_json) => {
        Plotly.newPlot('analysePlot', response_json);
        document.getElementById('analysePlot').on('plotly_click', (data) => {
          console.log(data);
          config.set_step(data.points[0].pointIndex);
        });
      });
  };
}

function clickAddSceneModifier() {
  document.getElementById('addSceneModifierImportBtn').onclick = async function () {
    const function_id = document.getElementById(
      'addSceneModifierImport',
    ).value;
    await addSceneModifierOption(function_id);
    document
      .getElementById('addSceneModifier')
      .dispatchEvent(new Event('change'));
  };
}

function resizeOffcanvas() {
  // Rescale offcanvas by dragging
  let active_offcanvas_border;
  const offcanvas_borders = document.getElementsByClassName('offcanvas-border');

  function resize_offcanvas(e) {
    if (e.clientX < 200) {
      return;
    }
    active_offcanvas_border.parentNode.style.width = `${e.clientX}px`;
    active_offcanvas_border.style.left = `${e.clientX}px`;
  }

  for (let i = 0; i < offcanvas_borders.length; i++) {
    offcanvas_borders[i].onpointerdown = function (e) {
      console.log(this);
      active_offcanvas_border = this;
      document.addEventListener('pointermove', resize_offcanvas);
    };
  }

  document.addEventListener('pointerup', (e) => {
    document.removeEventListener('pointermove', resize_offcanvas);
  });
}

export function setUIEvents(socket, world) {
  // setupSlider(socket, world);
  // setupPlayPause(world);
  // update_materials(config);
  // updateFPS(config);
  // setupSlider(config);
  // setupPlayPause(config);
  // attachKeyPressed(config);
  // update_resolution(config, world);
  // update_sphere_radius(config);
  // update_bond_radius(config, world);
  // loadSceneModifier(config, world);
  // loadSceneAnalysis(config, world);

  // clickAddSceneModifier();
  resizeOffcanvas();

  // // disable loading spinner by making it invisible
  // const loadingElem = document.getElementById("atom-spinner");
  // loadingElem.style.display = "none";
}
