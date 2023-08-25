JSONEditor.defaults.options.theme = "bootstrap5";
JSONEditor.defaults.options.iconlib = "fontawesome5";
JSONEditor.defaults.options.object_background = "bg-white";
JSONEditor.defaults.options.disable_edit_json = true;
JSONEditor.defaults.options.disable_properties = true;
JSONEditor.defaults.options.disable_collapse = true;

export function initJSONEditor(socket, cache, world) {
  let editor = new JSONEditor(
    document.getElementById("interaction-json-editor"),
    {
      schema: { type: "object", title: "ZnDraw", properties: {} },
    },
  );
  socket.on("modifier:schema", function (data) {
    editor.destroy();
    editor = new JSONEditor(
      document.getElementById("interaction-json-editor"),
      {
        schema: data,
      },
    );
  });
  document
    .getElementById("interaction-json-editor-submit")
    .addEventListener("click", function () {
      // Get the value from the editor
      const value = editor.getValue();
      console.log(value);
      socket.emit("modifier:run", {
        "params": value,
        "atoms": cache.get(world.getStep()),
        "selection": world.getSelection(),
        "step": world.getStep(),
      });
    });
  socket.emit("modifier:schema");
}
