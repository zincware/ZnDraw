JSONEditor.defaults.options.theme = 'bootstrap5';
JSONEditor.defaults.options.iconlib = 'fontawesome5';
JSONEditor.defaults.options.object_background = 'bg-white';
JSONEditor.defaults.options.disable_edit_json = true;
JSONEditor.defaults.options.disable_properties = true;
JSONEditor.defaults.options.disable_collapse = true;
JSONEditor.defaults.options.no_additional_properties = true;

export function initJSONEditor(socket, cache, world) {
  modifier_editor(socket, cache, world);
  analysis_editor(socket, cache, world);
  scene_editor(socket, cache, world);
  selection_editor(socket, cache, world);
}

function selection_editor(socket, cache, world) {
  const selection = document.getElementById('selection-select');

  socket.on('selection:schema', (data) => {
    const option = document.createElement('option');
    option.value = JSON.stringify(data.schema);
    option.innerHTML = data.name;
    selection.appendChild(option);
  });

  socket.emit('selection:schema');

  let editor;
  selection.onchange = function () {
    if (editor !== undefined) {
      editor.destroy();
    }
    if (selection.value === '') {
      return;
    }
    const schema = JSON.parse(selection.value);
    editor = new JSONEditor(
      document.getElementById('selection-json-editor'),
      {
        schema,
      },
    );
  };

  document.getElementById('selection-json-editor-submit').addEventListener('click', () => {
    // Get the value from the editor
    const value = editor.getValue();
    console.log(value);

    socket.emit('selection:run', {
      name: selection.options[selection.selectedIndex].text,
      params: value,
      atoms: cache.get(world.getStep()),
      selection: world.getSelection(),
    });
  });
}

function scene_editor(socket, cache, world) {
  socket.emit('scene:schema', (data) => {
    const editor = new JSONEditor(
      document.getElementById('scene-json-editor'),
      {
        schema: data,
      },
    );
    editor.on('change', () => {
      const value = editor.getValue();
      world.rebuild(value.resolution, value.material, value.wireframe, value.simulation_box, value.bonds);
      world.player.setLoop(value['Animation Loop']);
    });
  });
}

function analysis_editor(socket, cache, world) {
  let editor;
  const selection = document.getElementById('analysis-select');
  selection.onchange = function () {
    if (editor !== undefined) {
      editor.destroy();
    }
    if (selection.value === '') {
      return;
    }
    const schema = JSON.parse(selection.value);
    editor = new JSONEditor(
      document.getElementById('analysis-json-editor'),
      {
        schema,
      },
    );
  };

  socket.on('analysis:schema', (data) => {
    const option = document.createElement('option');
    option.value = JSON.stringify(data.schema);
    option.innerHTML = data.name;
    selection.appendChild(option);
  });

  document.getElementById('analysis-json-editor-submit').addEventListener('click', () => {
    // Get the value from the editor
    const value = editor.getValue();
    console.log(value);

    socket.emit('analysis:run', {
      name: selection.options[selection.selectedIndex].text,
      params: value,
      atoms: cache.get(world.getStep()),
      selection: world.getSelection(),
      step: world.getStep(),
      atoms_list: cache.getAllAtoms(),
    }, (data) => {
      Plotly.newPlot('analysisPlot', JSON.parse(data));

      function buildPlot() {
        Plotly.newPlot('analysisPlot', JSON.parse(data));
        const myplot = document.getElementById('analysisPlot');
        myplot.on('plotly_click', (data) => {
          const point = data.points[0];
          const step = point.x;
          world.setStep(step);
        });
      }

      // resizeObserver = new ResizeObserver(() => {
      //   buildPlot();
      // }).observe(document.getElementById("analysisPlot").parentElement);

      buildPlot();
    });

    document.getElementById('analysis-json-editor-submit').disabled = true;
    // if there is an error in uploading, we still want to be able to submit again
    setTimeout(() => {
      document.getElementById('analysis-json-editor-submit').disabled = false;
    }, 1000);
  });

  function get_analysis_data() {
    if (cache.get(0) !== undefined) {
      socket.emit('analysis:schema', { atoms: cache.get(0) });
    } else {
      setTimeout(get_analysis_data, 100);
    }
  }
  get_analysis_data();
}

function modifier_editor(socket, cache, world) {
  let editor;
  const selection = document.getElementById('modifier-select');

  selection.onchange = function () {
    if (editor !== undefined) {
      editor.destroy();
    }
    if (selection.value === '') {
      return;
    }
    const schema = JSON.parse(selection.value);
    editor = new JSONEditor(
      document.getElementById('interaction-json-editor'),
      {
        schema,
      },
    );
  };

  socket.on('modifier:schema', (data) => {
    const option = document.createElement('option');
    option.value = JSON.stringify(data.schema);
    option.innerHTML = data.name;
    selection.appendChild(option);
  });

  document
    .getElementById('interaction-json-editor-submit')
    .addEventListener('click', () => {
      // Get the value from the editor
      const value = editor.getValue();
      console.log(value);

      const { points, segments } = world.getLineData();

      socket.emit('modifier:run', {
        name: selection.options[selection.selectedIndex].text,
        params: value,
        atoms: cache.get(world.getStep()),
        selection: world.getSelection(),
        step: world.getStep(),
        points,
        segments,
      });

      document.getElementById('interaction-json-editor-submit').disabled = true;
      // if there is an error in uploading, we still want to be able to submit again
      setTimeout(() => {
        document.getElementById('interaction-json-editor-submit').disabled = false;
      }, 1000);
    });

  socket.emit('modifier:schema');
}
