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

export function duplicateAnchorPoints() {
  // TODO shift added point a bit + insert at correct position
  const transform_object = this.transform_controls.object;
  if (transform_object.name === "AnchorPoint" && !this._drawing) {
    const index = this.line3D.anchorPoints.children.indexOf(transform_object);

    let point;
    if (index > 0) {
      // Add the point between the current and the previous point
      const obj_before = this.line3D.anchorPoints.children[index - 1];
      const new_pos = obj_before.position
        .clone()
        .sub(transform_object.position)
        .multiplyScalar(0.5)
        .add(transform_object.position);

      point = this.line3D.addPoint(new_pos, index);
    } else {
      if (this.line3D.anchorPoints.children.length > 1) {
        // Add the point between the current and the next point
        const obj_after = this.line3D.anchorPoints.children[index + 1];
        const new_pos = obj_after.position
          .clone()
          .sub(transform_object.position)
          .multiplyScalar(0.5)
          .add(transform_object.position);

        point = this.line3D.addPoint(new_pos, index - 1);
      }
    }
    this.transform_controls.detach();
    this.transform_controls.attach(point);
  }
}
