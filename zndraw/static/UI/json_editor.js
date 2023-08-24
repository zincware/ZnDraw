export function initJSONEditor() {
    JSONEditor.defaults.options.theme = 'bootstrap5';
    JSONEditor.defaults.options.iconlib = 'fontawesome5';
    JSONEditor.defaults.options.object_background = 'bg-white';
    JSONEditor.defaults.options.disable_edit_json = true;
    JSONEditor.defaults.options.disable_properties = true;
    JSONEditor.defaults.options.disable_collapse = true;

    // Initialize the editor with a JSON schema
    var editor = new JSONEditor(document.getElementById('interaction-json-editor'), 
    {
        schema: {
            type: "object",
            title: "Car",
            properties: {
                safety: {
                    type: "integer",
                    minimum: "0",
                    maximum: "5",
                    format: "range",
                }
            }
        }
    });

    // Hook up the submit button to log to the console
    document.getElementById('interaction-json-editor-submit').addEventListener('click', function () {
        // Get the value from the editor
        console.log(editor.getValue());
    });
}