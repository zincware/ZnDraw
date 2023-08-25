const setSize = (container, camera, renderer) => {
  camera.aspect = container.clientWidth / container.clientHeight;
  camera.updateProjectionMatrix();

  renderer.setSize(container.clientWidth, container.clientHeight);
  renderer.setPixelRatio(window.devicePixelRatio);
};

class Resizer {
  constructor(container, camera, renderer, renderer2d) {
    // set initial size on load
    setSize(container, camera, renderer);

    window.addEventListener("resize", () => {
      // set the size again if a resize occurs
      setSize(container, camera, renderer);
      renderer2d.setSize(window.innerWidth, window.innerHeight);
      // perform any custom actions
      this.onResize();
    });
  }

  onResize() {}
}

export { Resizer };
