// Interaction with Python methods
import { createElementFromSchema } from './schemaforms.js';
import * as DATA from './data.js';

export const addModifierModal = new bootstrap.Modal(document.getElementById("addModifierModal"));
export const addAnalysisModal = new bootstrap.Modal(document.getElementById("addAnalysisModal"));

export async function addAnalysisOption(function_id) {
    await fetch("add_analysis", {
        "method": "POST",
        "headers": { "Content-Type": "application/json" },
        "body": JSON.stringify(function_id),
    }).then(response => response.json()).then(function (response_json) {
        // if not null alert
        if ("error" in response_json) {
            // TODO check if method is already loaded
            alert(response_json["error"]);
            stepError(response_json["error"]);
        } else {
            if (document.getElementById("scene-analysis_" + response_json["title"]) != null) {
                alert("Function already loaded");
                stepError("Function already loaded");
            };
            addAnalysisModal.hide();
            DATA.load_config();
        }
        return response_json;
    }).then(function (response_json) {
        let modifier = document.createElement("option");
        modifier.value = response_json["title"];
        modifier.innerHTML = response_json["title"];
        document.getElementById("addAnalysis").appendChild(modifier);
        return response_json;
    }).then(function (response_json) {
        let analysisSettings = document.getElementById("analysisSettings");
        analysisSettings.appendChild(createElementFromSchema(response_json, "scene-analysis"));
        document.getElementById("addAnalysis").value = response_json["title"];
    });
}

export async function addSceneModifierOption(function_id) {
    await fetch("add_update_function", {
        "method": "POST",
        "headers": { "Content-Type": "application/json" },
        "body": JSON.stringify(function_id),
    }).then(response => response.json()).then(function (response_json) {
        // if not null alert
        if ("error" in response_json) {
            // TODO check if method is already loaded
            alert(response_json["error"]);
            stepError(response_json["error"]);
        } else {
            if (document.getElementById("scene-modifier_" + response_json["title"]) != null) {
                alert("Function already loaded");
                stepError("Function already loaded");
            };
            addModifierModal.hide();
            DATA.load_config();
        }
        return response_json;
    }).then(function (response_json) {
        let modifier = document.createElement("option");
        modifier.value = response_json["title"];
        modifier.innerHTML = response_json["title"];
        document.getElementById("addSceneModifier").appendChild(modifier);
        return response_json;
    }).then(function (response_json) {
        let sceneModifierSettings = document.getElementById("sceneModifierSettings");
        sceneModifierSettings.appendChild(createElementFromSchema(response_json, "scene-modifier"));
        document.getElementById("addSceneModifier").value = response_json["title"];
    });
}


// load analysis methods from config 
export async function loadAnalysisMethods() {
    // iterate DATA.config.analysis_methods and add them to the select
    for (let i = 0; i < DATA.config.analysis_functions.length; i++) {
        let analysis_function = DATA.config.analysis_functions[i];
        await addAnalysisOption(analysis_function);
    }
    document.getElementById("addAnalysis").value = "";
    document.getElementById("addAnalysis").dispatchEvent(new Event('change'));
}

// load analysis methods from config 
export async function loadSceneModifier() {
    // iterate DATA.config.analysis_methods and add them to the select
    for (let i = 0; i < DATA.config.modify_functions.length; i++) {
        let modify_function = DATA.config.modify_functions[i];
        await addSceneModifierOption(modify_function);
    }
    document.getElementById("addSceneModifier").value = "";
    document.getElementById("addSceneModifier").dispatchEvent(new Event('change'));
}
