JSONEditor.defaults.options.theme = "bootstrap5";
JSONEditor.defaults.options.iconlib = "fontawesome5";
JSONEditor.defaults.options.object_background = "bg-body-color";
JSONEditor.defaults.options.disable_edit_json = true;
JSONEditor.defaults.options.disable_properties = true;
JSONEditor.defaults.options.disable_collapse = true;
JSONEditor.defaults.options.no_additional_properties = true;
JSONEditor.defaults.options.keep_oneof_values = false;

export function initJSONEditor(socket, cache, world) {
  modifier_editor(socket, cache, world);
  analysis_editor(socket, cache, world);
  scene_editor(socket, cache, world);
  selection_editor(socket, cache, world);
  draw_editor(socket, cache, world);
}

function rebuildEditor(editor, localStorageKey, data, div) {
  if (editor) {
    editor.destroy();
  }
  editor = new JSONEditor(div, {
    schema: data,
  });
  editor.on("change", () => {
    const editorValue = editor.getValue();
    editor.validate();
    localStorage.setItem(localStorageKey, JSON.stringify(editorValue));
  });

  editor.on("ready", () => {
    const userInput = localStorage.getItem(localStorageKey);
    if (userInput) {
      editor.setValue(JSON.parse(userInput));
    }
  });
  return editor;
}

function setupBtnQueue(socket, btn, task) {
  // document.getElementById("analysis-json-editor-submit")
  const originalInnerHTML = btn.innerHTML;

  socket.on(`${task}:run:running`, () => {
    btn.innerHTML = '<i class="fa-solid fa-spinner"></i> Running';
  });

  // Finished running
  socket.on(`${task}:run:finished`, (data) => {
    btn.innerHTML = originalInnerHTML;
    btn.disabled = false;
  });

  // other client started process
  socket.on(`${task}:run:enqueue`, () => {
    btn.disabled = true;
    btn.innerHTML = '<i class="fa-solid fa-hourglass-start"></i> Job queued';
  });

  // TODO: queue update
  // socket.on("modifier:queue:update", (position) => {
  //   document.getElementById("interaction-json-editor-submit").disabled = true;
  //   document.getElementById("interaction-json-editor-submit").innerHTML =
  //     '<i class="fa-solid fa-hourglass-start"></i> Job queued at position ' +
  //     position;
  //   // emit event "modifier:queue:update"
  //   const event = new CustomEvent("modifier:queue:update", {
  //     detail: { position: position },
  //   });
  //   document.dispatchEvent(event);
  // });
}

function draw_editor(socket, cache, world) {
  let editor;
  const btn = document.getElementById("drawCardBtn");
  const div = document.getElementById("draw-json-editor");

  btn.addEventListener("click", () => {
    socket.emit("draw:schema", (data) => {
      editor = rebuildEditor(editor, "draw-json-editor-input", data, div);
      // TODO: the button handling is all over the place - this should also be done through / with Python
      editor.on("change", () => {
        document.getElementById("drawAddCanvas").parameters = editor.getValue();
      });
    });
  });
}

function selection_editor(socket, cache, world) {
  let editor;
  const btn = document.getElementById("selectionCardBtn");
  const div = document.getElementById("selection-json-editor");

  btn.addEventListener("click", () => {
    socket.emit("selection:schema", (data) => {
      editor = rebuildEditor(editor, "selection-json-editor-input", data, div);
      // TODO: the button handling is all over the place - this should also be done through / with Python
      editor.on("change", () => {
        const value = editor.getValue();
        div.parameters = value;
      });
    });
  });

  setupBtnQueue(
    socket,
    document.getElementById("selection-json-editor-submit"),
    "selection",
  );

  document
    .getElementById("selection-json-editor-submit")
    .addEventListener("click", () => {
      // console.log(new Date().toISOString(), "running selection");
      // Get the value from the editor
      const errors = editor.validate();
      if (errors.length) {
        console.log(errors);
      } else {
        const value = editor.getValue();
        socket.emit("selection:run", value);
      }
    });
}

function scene_editor(socket, cache, world) {
  const btn = document.getElementById("sceneCardBtn");
  const div = document.getElementById("scene-json-editor");
  let editor;

  btn.addEventListener("click", () => {
    socket.emit("scene:schema", (data) => {
      editor = rebuildEditor(editor, "scene-json-editor-input", data, div);
      // TODO: the button handling is all over the place - this should also be done through / with Python
      editor.on("change", () => {
        const value = editor.getValue();
        world.rebuild(
          value.resolution,
          value.material,
          value.wireframe,
          value.simulation_box,
          value.bonds,
          value.label_offset,
          value.particle_size,
          value.bonds_size,
          value.fps,
          value.line_label,
        );
        world.player.setLoop(value["Animation Loop"]);
      });
    });
  });
}

function analysis_editor(socket, cache, world) {
  const div = document.getElementById("analysis-json-editor");
  const btn = document.getElementById("analysisCardBtn");
  let editor;

  btn.addEventListener("click", () => {
    socket.emit("analysis:schema", (data) => {
      editor = rebuildEditor(editor, "analysis-json-editor-input", data, div);
    });
  });

  setupBtnQueue(
    socket,
    document.getElementById("analysis-json-editor-submit"),
    "analysis",
  );

  document
    .getElementById("analysis-json-editor-submit")
    .addEventListener("click", () => {
      // Get the value from the editor
      const errors = editor.validate();
      if (errors.length) {
        console.log(errors);
      } else {
        const value = editor.getValue();
        // responseReceived = false;

        socket.emit("analysis:run", value);

        document.getElementById("analysis-json-editor-submit").disabled = true;
      }
    });

  socket.on("analysis:figure:set", (data) => {
    console.log(data);
    Plotly.newPlot("analysisPlot", JSON.parse(data));

    function buildPlot() {
      Plotly.newPlot("analysisPlot", JSON.parse(data));
      const myplot = document.getElementById("analysisPlot");
      myplot.on("plotly_click", (data) => {
        const point = data.points[0];
        const step = point.x;
        world.setStep(step);
      });
    }

    buildPlot();
  });
}

function modifier_editor(socket, cache, world) {
  const div = document.getElementById("interaction-json-editor");
  const btn = document.getElementById("interactionCardBtn");

  let editor = new JSONEditor(div, {
    schema: { type: "object", title: "ZnDraw", properties: {} },
  });

  socket.on("modifier:queue:update", (position) => {
    document.getElementById("interaction-json-editor-submit").disabled = true;
    document.getElementById("interaction-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-hourglass-start"></i> Job queued at position ' +
      position;
    // emit event "modifier:queue:update"
    const event = new CustomEvent("modifier:queue:update", {
      detail: { position: position },
    });
    document.dispatchEvent(event);
  });

  setupBtnQueue(
    socket,
    document.getElementById("interaction-json-editor-submit"),
    "modifier",
  );

  // // Check if a running response is received
  // socket.on("modifier:run:running", () => {
  //   document.getElementById("interaction-json-editor-submit").innerHTML =
  //     '<i class="fa-solid fa-spinner"></i> Running';

  //   const event = new CustomEvent("modifier:run:running");
  //   document.dispatchEvent(event);
  // });

  // // Finished running
  // socket.on("modifier:run:finished", () => {
  //   document.getElementById("interaction-json-editor-submit").innerHTML =
  //     '<i class="fa-solid fa-play"></i> Run Modifier';
  //   document.getElementById("interaction-json-editor-submit").disabled = false;
  //   const event = new CustomEvent("modifier:run:finished");
  //   document.dispatchEvent(event);
  // });

  // // other client started process
  // socket.on("modifier:run:enqueue", () => {
  //   document.getElementById("interaction-json-editor-submit").disabled = true;
  //   document.getElementById("interaction-json-editor-submit").innerHTML =
  //     '<i class="fa-solid fa-hourglass-start"></i> Job queued';
  //   const event = new CustomEvent("modifier:run:enqueue");
  //   document.dispatchEvent(event);
  // });

  document
    .getElementById("interaction-json-editor-submit")
    .addEventListener("click", () => {
      // Get the value from the editor
      const errors = editor.validate();
      if (errors.length) {
        console.log(errors);
      } else {
        const value = editor.getValue();
        // responseReceived = false;

        socket.emit("modifier:run", value);
        const event = new CustomEvent("modifier:run");
        document.dispatchEvent(event);

        document.getElementById("interaction-json-editor-submit").disabled =
          true;
      }
    });

  // create when btn is pressed
  btn.addEventListener("click", () => {
    socket.emit("modifier:schema", (data) => {
      editor = rebuildEditor(
        editor,
        "interaction-json-editor-input",
        data,
        div,
      );
    });
  });
}
