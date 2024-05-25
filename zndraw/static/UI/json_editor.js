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

function draw_editor(socket, cache, world) {
  socket.on("draw:schema", (data) => {
    const div = document.getElementById("draw-json-editor");
    const editor = new JSONEditor(div, {
      schema: data,
    });
    editor.on("change", () => {
      document.getElementById("drawAddCanvas").parameters = editor.getValue();
    });
  });
}

function selection_editor(socket, cache, world) {
  socket.on("selection:schema", (data) => {
    const div = document.getElementById("selection-json-editor");
    const editor = new JSONEditor(div, {
      schema: data,
    });

    editor.on("change", () => {
      const value = editor.getValue();
      div.parameters = value;
    });

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
  });
}

function scene_editor(socket, cache, world) {
  socket.on("scene:schema", (data) => {
    const editor = new JSONEditor(
      document.getElementById("scene-json-editor"),
      {
        schema: data,
      },
    );
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
}

function analysis_editor(socket, cache, world) {
  const div = document.getElementById("analysis-json-editor");
  let editor = new JSONEditor(div, {
    schema: { type: "object", title: "Analysis", properties: {} },
  });

  socket.on("analysis:schema", (data) => {
    editor.destroy();
    editor = new JSONEditor(div, {
      schema: data,
    });
  });

  socket.on("analysis:run:running", () => {
    document.getElementById("analysis-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-spinner"></i> Running';
  });

  // Finished running
  socket.on("analysis:run:finished", (data) => {
    document.getElementById("analysis-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-play"></i> Analyse';
    document.getElementById("analysis-json-editor-submit").disabled = false;
  });

  // other client started process
  socket.on("analysis:run:enqueue", () => {
    document.getElementById("analysis-json-editor-submit").disabled = true;
    document.getElementById("analysis-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-hourglass-start"></i> Job queued';
  });

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

  socket.on("analysis:figure", (data) => {
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

  // Check if a running response is received
  socket.on("modifier:run:running", () => {
    document.getElementById("interaction-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-spinner"></i> Running';

    const event = new CustomEvent("modifier:run:running");
    document.dispatchEvent(event);
  });

  // Finished running
  socket.on("modifier:run:finished", () => {
    document.getElementById("interaction-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-play"></i> Run Modifier';
    document.getElementById("interaction-json-editor-submit").disabled = false;
    const event = new CustomEvent("modifier:run:finished");
    document.dispatchEvent(event);
  });

  // other client started process
  socket.on("modifier:run:enqueue", () => {
    document.getElementById("interaction-json-editor-submit").disabled = true;
    document.getElementById("interaction-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-hourglass-start"></i> Job queued';
    const event = new CustomEvent("modifier:run:enqueue");
    document.dispatchEvent(event);
  });

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
      const localStorageKey = "interaction-json-editor-input"
      
      editor.destroy();
      editor = new JSONEditor(div, {
        schema: data,
      });
      editor.on('change', () => {
        const editorValue = editor.getValue();
        localStorage.setItem(localStorageKey, JSON.stringify(editorValue));
      });

      editor.on('ready',() => {
        const userInput = localStorage.getItem(localStorageKey);
        if (userInput) {
          editor.setValue(JSON.parse(userInput));
        };        
      });
      
    });
  });
}
