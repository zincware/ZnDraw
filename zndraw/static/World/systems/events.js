import * as THREE from "three";

export function centerCamera(controls, particlesGroup) {
  if (controls.enablePan) {
    // get the first object that is selected
    const dummy = new THREE.Vector3();
    dummy.copy(controls.getCenter(particlesGroup.selection));
    controls.target.copy(dummy);
  } else {
    // follow is currently not working due to instancing
    // document.getElementById("alertBoxCamera").style.display = "none";
    // controls.target = controls.target.clone();
    // controls.enablePan = true;
  }
}
