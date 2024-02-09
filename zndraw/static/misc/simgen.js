const setupSiMGen = function (socket, world) {
  const runButton = document.getElementById("simgen-run-btn");
  const linkButton = document.getElementById("simgen-link-btn");
  const canvasButton = document.getElementById("simgen-canvas-btn");
  if (!runButton || !linkButton || !canvasButton) {
    console.log("simgen buttons not found");
    return;
  }

  const runButtonText = runButton.innerHTML;
  const linkButtonText = linkButton.innerHTML;
  const canvasButtonText = canvasButton.innerHTML;

  let clickedButton;

  runButton.addEventListener("click", () => {
    clickedButton = runButton;

    runButton.disabled = true;
    linkButton.disabled = true;
    canvasButton.disabled = true;
    socket.emit("modifier:run", {
      method: { discriminator: "Delete" },
    });
  });

  linkButton.addEventListener("click", () => {
    clickedButton = linkButton;

    runButton.disabled = true;
    linkButton.disabled = true;
    canvasButton.disabled = true;
    socket.emit("modifier:run", {
      method: { discriminator: "Delete" },
    });
  });

  canvasButton.addEventListener("click", () => {
    // clickedButton = canvasButton;
    // runButton.disabled = true;
    // linkButton.disabled = true;
    // canvasButton.disabled = true;
    // TODO: shouldn't trash just be another modifier?
    socket.emit("scene:trash");
    // clickt drawAddCanvas btn
    document.getElementById("drawAddCanvas").click();
    // click after some time, trash needs to be done first
    setTimeout(function () {
      document.getElementById("drawingSwitch").click();
    }, 1000);
  });

  // act on custom event "modifier:queue:update"
  document.addEventListener("modifier:queue:update", (event) => {
    clickedButton.innerHTML =
      '<i class="fa-solid fa-hourglass-start"></i> Job queued at position ' +
      event.detail.position;
  });

  document.addEventListener("modifier:run:enqueue", () => {
    clickedButton.innerHTML =
      '<i class="fa-solid fa-hourglass-start"></i> Job queued';
  });

  document.addEventListener("modifier:run:finished", () => {
    if (clickedButton === runButton) {
      runButton.innerHTML = runButtonText;
    } else if (clickedButton === linkButton) {
      linkButton.innerHTML = linkButtonText;
    } else if (clickedButton === canvasButton) {
      canvasButton.innerHTML = canvasButtonText;
    }
    runButton.disabled = false;
    linkButton.disabled = false;
    canvasButton.disabled = false;
    clickedButton = null;
  });
};

export { setupSiMGen };
