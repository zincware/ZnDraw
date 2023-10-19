JSONEditor.defaults.options.theme = "bootstrap5";
JSONEditor.defaults.options.iconlib = "fontawesome5";
JSONEditor.defaults.options.object_background = "bg-white";
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
        console.log(new Date().toISOString(), "running selection");
        // Get the value from the editor
        const value = editor.getValue();
        console.log(value);

        socket.emit("selection:run", {
          params: value,
          atoms: cache.get(world.getStep()),
          selection: world.getSelection(),
        });
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

        document.getElementById("analysis-json-editor-submit").disabled = true;
        // if there is an error in uploading, we still want to be able to submit again
        setTimeout(() => {
          document.getElementById(
            "analysis-json-editor-submit",
          ).disabled = false;
        }, 1000);
      });
  });
}

function modifier_editor(socket, cache, world) {
  socket.on("modifier:schema", (data) => {
    const div = document.getElementById("interaction-json-editor");
    const editor = new JSONEditor(div, {
      schema: data,
    });

    editor.on("change", () => {
      const value = editor.getValue();
      div.parameters = value;
    });

    document
      .getElementById("interaction-json-editor-submit")
      .addEventListener("click", () => {
        // Get the value from the editor
        const value = editor.getValue();
        const { points, segments } = world.getLineData();
        console.log(value);
        socket.emit("modifier:run", {
          params: value,
          url: window.location.href,
        });
        // world.particles.click(); // reset selection

        document.getElementById(
          "interaction-json-editor-submit",
        ).disabled = true;
        // if there is an error in uploading, we still want to be able to submit again
        setTimeout(() => {
          document.getElementById(
            "interaction-json-editor-submit",
          ).disabled = false;
        }, 1000);
      });
  });
}
