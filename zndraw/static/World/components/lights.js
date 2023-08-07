import * as THREE from "three";

function createLights() {
  // Create a directional light
  const hemisphereLight = new THREE.HemisphereLight(0xffffff, 0x777777, 0.1);
  const spotLight = new THREE.SpotLight(0xffffff, 0, 0, Math.PI / 2);
  return [hemisphereLight, spotLight];
}

export { createLights };
