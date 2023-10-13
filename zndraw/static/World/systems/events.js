import * as THREE from "three";

export function centerCamera(controls, particlesGroup) {
  if (controls.enablePan) {
    // get the first object that is selected
    if (particlesGroup.selection.length > 0) {
      const matrix = new THREE.Matrix4();
      const dummy = new THREE.Object3D();
      particlesGroup.particles_mesh.getMatrixAt(
        particlesGroup.selection[0],
        matrix,
      );
      matrix.decompose(dummy.position, dummy.quaternion, dummy.scale);
      controls.target.copy(dummy.position);
      // controls.enablePan = false;
      // document.getElementById("alertBoxCamera").style.display = "block";
    } else {
      const dummy = new THREE.Vector3();
      dummy.copy(controls.getCenter());
      controls.target.copy(dummy);
    }
  } else {
    // follow is currently not working due to instancing
    // document.getElementById("alertBoxCamera").style.display = "none";
    // controls.target = controls.target.clone();
    // controls.enablePan = true;
  }
}
