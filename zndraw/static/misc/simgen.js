const setupSiMGen = function (socket, world) {
  const runButton = document.getElementById("simgen-run-btn");
  const linkButton = document.getElementById("simgen-link-btn");
  const canvasButton = document.getElementById("simgen-canvas-btn");
  if (!runButton || !linkButton || !canvasButton) {
    console.log("simgen buttons not found");
    return;
  }
  const consoleBtn = document.getElementById("ZnDrawConsoleSwitch");
  consoleBtn.click();

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
      method: { discriminator: "SiMGenDemo" },
    });
  });

  linkButton.addEventListener("click", () => {
    clickedButton = linkButton;

    runButton.disabled = true;
    linkButton.disabled = true;
    canvasButton.disabled = true;
    socket.emit("modifier:run", {
      method: { discriminator: "Connect" },
    });
  });

  canvasButton.addEventListener("click", () => {
    socket.emit("modifier:run", {
      method: { discriminator: "ClearScene" },
    });
    setTimeout(function () {
      document.getElementById("drawAddCanvas").click();
      document.getElementById("drawingSwitch").click();
    }, 1000);
  });

  document.addEventListener("modifier:run:running", () => {
    clickedButton.innerHTML = '<i class="fa-solid fa-spinner"></i> Running';
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

  document.addEventListener("modifier:run", () => {
    runButton.disabled = false;
    linkButton.disabled = false;
    canvasButton.disabled = false;
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
