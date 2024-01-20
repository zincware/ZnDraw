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
          console.log(value);

          socket.emit("selection:run", {
            params: value,
          });
        }
      });
  });
}

function scene_editor(socket, cache, world) {
  socket.emit("scene:schema", (data) => {
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
  socket.on("analysis:schema", (data) => {
    const div = document.getElementById("analysis-json-editor");
    const editor = new JSONEditor(div, {
      schema: data,
    });

    editor.on("change", () => {
      const value = editor.getValue();
      div.parameters = value;
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

          socket.emit("analysis:run", {
            params: value,
          });

          document.getElementById("analysis-json-editor-submit").disabled =
            true;
          // if there is an error in uploading, we still want to be able to submit again
          setTimeout(() => {
            document.getElementById("analysis-json-editor-submit").disabled =
              false;
          }, 1000);
        }
      });
  });
}

function modifier_editor(socket, cache, world) {
  const div = document.getElementById("interaction-json-editor");
  let editor = new JSONEditor(div, {
    schema: { type: "object", title: "ZnDraw", properties: {} },
  });

  let responseReceived = false;

  socket.on("modifier:run:submitted", () => {
    console.log(new Date().toISOString(), "modifier:run:submitted");
    setTimeout(() => {
      if (!responseReceived) {
        console.warn("No response on 'modifier:run:running' within 1 second");
        alert("No response from server. Please try again.");
        document.getElementById("interaction-json-editor-submit").disabled =
          false;
        socket.emit("modifier:run:failed");
      }
    }, 1000);
  });

  socket.on("modifier:run:enqueue", (position) => {
    console.log(new Date().toISOString(), "modifier:run:enqueue");
    document.getElementById("interaction-json-editor-submit").disabled = true;
    document.getElementById("interaction-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-hourglass-start"></i> Job queued at position ' +
      position;
  });

  // Check if a running response is received
  socket.on("modifier:run:running", () => {
    responseReceived = true;
    document.getElementById("interaction-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-spinner"></i> Running';
  });

  // Finished running
  socket.on("modifier:run:finished", (data) => {
    console.log(new Date().toISOString(), "modifier:run:finished");
    document.getElementById("interaction-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-play"></i> Run Modifier';
    document.getElementById("interaction-json-editor-submit").disabled = false;
  });

  socket.on("modifier:run:criticalfail", (data) => {
    console.log(new Date().toISOString(), "modifier:run:failed");
    alert("Modifier failed. Please try again with different settings.");
    document.getElementById("interaction-json-editor-submit").innerHTML =
      '<i class="fa-solid fa-play"></i> Run Modifier';
    document.getElementById("interaction-json-editor-submit").disabled = false;
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
        console.log(value);

        responseReceived = false;

        socket.emit("modifier:run", {
          params: value,
          url: window.location.href,
        });

        document.getElementById("interaction-json-editor-submit").disabled =
          true;
      }
    });

  socket.on("modifier:schema", (data) => {
    editor.destroy();

    editor = new JSONEditor(div, {
      schema: data,
    });

    editor.on("change", () => {
      const value = editor.getValue();
      div.parameters = value;
    });
  });
}
