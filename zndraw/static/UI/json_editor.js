JSONEditor.defaults.options.theme = 'bootstrap5';
JSONEditor.defaults.options.iconlib = 'fontawesome5';
JSONEditor.defaults.options.object_background = 'bg-white';
JSONEditor.defaults.options.disable_edit_json = true;
JSONEditor.defaults.options.disable_properties = true;
JSONEditor.defaults.options.disable_collapse = true;
JSONEditor.defaults.options.no_additional_properties = true;

export function initJSONEditor(socket, cache, world) {
  let editor = new JSONEditor(
    document.getElementById('interaction-json-editor'),
    {
      schema: { type: 'object', title: 'ZnDraw', properties: {} },
    },
  );
  const selection = document.getElementById('modifier-select');

  selection.onchange = function () {
    const schema = JSON.parse(selection.value);
    editor.destroy();
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
    // select the one that was just added
    selection.value = JSON.stringify(data.schema);

    editor.destroy();
    editor = new JSONEditor(
      document.getElementById('interaction-json-editor'),
      {
        schema: data.schema,
      },
    );
  });

  document
    .getElementById('interaction-json-editor-submit')
    .addEventListener('click', () => {
      // Get the value from the editor
      const value = editor.getValue();
      console.log(value);

      const {points, segments}  = world.getLineData();

      socket.emit('modifier:run', {
        name: selection.options[selection.selectedIndex].text,
        params: value,
        atoms: cache.get(world.getStep()),
        selection: world.getSelection(),
        step: world.getStep(),
        points: points,
        segments: segments,
      });

      document.getElementById("interaction-json-editor-submit").disabled = true;
      // if there is an error in uploading, we still want to be able to submit again
      setTimeout(() => {
        document.getElementById("interaction-json-editor-submit").disabled = false;
      }, 1000);
    });

  socket.emit('modifier:schema');
}
