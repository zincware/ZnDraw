import { WebGLRenderer } from "three";
import { CSS2DRenderer } from "three/examples/jsm/renderers/CSS2DRenderer.js";

function createRenderer() {
  const renderer = new WebGLRenderer({ antialias: true });

  return renderer;
}

function create2DRenderer() {
  const labelRenderer = new CSS2DRenderer();
  labelRenderer.setSize(window.innerWidth, window.innerHeight);
  labelRenderer.domElement.style.position = "absolute";
  labelRenderer.domElement.style.top = "0px";
  labelRenderer.domElement.style.pointerEvents = "none";
  document
    .getElementById("scene-container")
    .appendChild(labelRenderer.domElement);

  return labelRenderer;
}

export { createRenderer, create2DRenderer };
