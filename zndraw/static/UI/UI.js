function setupUpload(socket) {
  const file = {
    dom: document.getElementById("fileInput"),
    binary: null,
  };

  const reader = new FileReader();

  // Because FileReader is asynchronous, store its
  // result when it finishes reading the file
  reader.addEventListener("load", () => {
    file.binary = reader.result;
    socket.emit("upload", {
      content: reader.result,
      filename: file.dom.files[0].name,
    });
  });

  // At page load, if a file is already selected, read it.
  if (file.dom.files[0]) {
    reader.readAsBinaryString(file.dom.files[0]);
  }

  // If not, read the file once the user selects it.
  file.dom.addEventListener("change", () => {
    if (reader.readyState === FileReader.LOADING) {
      reader.abort();
    }

    reader.readAsBinaryString(file.dom.files[0]);
  });
}

function setupMobile() {
  // set display style of .playPauseCtrl to block if on mobile
  if (
    /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(
      navigator.userAgent,
    )
  ) {
    document.getElementById("playPauseCtrl").style.display = "block";
  }
}

function setupDragDrop(socket) {
  const scene = document.getElementById("scene-container");

  scene.addEventListener("dragover", (event) => {
    event.preventDefault();
    // show the overlay as long as the file is dragged over the scene
    const overlay = document.getElementById("overlay");
    overlay.style.display = "block";
  });

  scene.addEventListener("dragleave", (event) => {
    event.preventDefault();
    // hide the overlay when the file is dragged out of the scene
    const overlay = document.getElementById("overlay");
    overlay.style.display = "none";
  });

  scene.addEventListener("drop", (event) => {
    event.preventDefault();
    // hide the overlay when the file is dropped
    const overlay = document.getElementById("overlay");
    overlay.style.display = "none";

    // read the file
    const file = event.dataTransfer.files[0];
    console.log(event);
    if (!file) {
      console.error("No file was dropped");
      return;
    }

    const reader = new FileReader();
    reader.readAsBinaryString(file);

    // send the file to the server
    reader.addEventListener("load", () => {
      socket.emit("upload", {
        content: reader.result,
        filename: file.name,
      });
    });
  });
}

function setupTrashClick(socket) {
  document.getElementById("trashBtn").addEventListener("click", () => {
    socket.emit("scene:trash", {});
  });
}

function setupNavbarLeft() {
  const tooltipTriggerList = document.querySelectorAll(
    '[data-bs-toggle="tooltip"]',
  );
  const tooltipList = [...tooltipTriggerList].map(
    (tooltipTriggerEl) => new bootstrap.Tooltip(tooltipTriggerEl),
  );

  let menu = document.querySelectorAll("[name=leftMenuInput]");
  // for each one, remove the check if is was checked before click
  const clickState = {};
  menu.forEach((item) => {
    clickState[item.id] = false;
    item.addEventListener("click", () => {
      if (clickState[item.id]) {
        item.checked = false;
        clickState[item.id] = false;
      } else {
        menu.forEach((item) => {
          clickState[item.id] = false;
        });

        clickState[item.id] = true;
      }
    });
  });
}

function switchColorScheme(world) {
  const switchBtn = document.getElementById("colorModeSwitch");

  switchBtn.addEventListener("click", () => {
    const theme = document.documentElement.getAttribute("data-bs-theme");
    if (theme === "dark") {
      document.documentElement.setAttribute("data-bs-theme", "light");
      world.scene.background.set(0xffffff);
      switchBtn.innerHTML = '<i class="fa-solid fa-sun"></i>';
    } else {
      document.documentElement.setAttribute("data-bs-theme", "dark");
      world.scene.background.set(0x000000);
      switchBtn.innerHTML = '<i class="fa-solid fa-moon"></i>';
    }
  });
}

function setupPointerFrameChange(world) {
  const progress = document.getElementById("frameProgress");

  progress.addEventListener("pointerdown", (event) => {
    // get the relative position of the pointer from left to right in 0..1
    const relativePosition = event.offsetX / window.innerWidth;
    const step = Math.round(relativePosition * progress.ariaValueMax);
    world.setStep(step);
  });
}

export function setUIEvents(socket, cache, world) {
  // resizeOffcanvas();
  setupUpload(socket);
  setupNavbarLeft();
  setupMobile();
  setupDragDrop(socket);
  setupTrashClick(socket);
  switchColorScheme(world);
  setupPointerFrameChange(world);

  socket.on("download:response", (data) => {
    const blob = new Blob([data], { type: "text/csv" });
    const elem = window.document.createElement("a");
    elem.href = window.URL.createObjectURL(blob);
    elem.download = "trajectory.xyz";
    document.body.appendChild(elem);
    elem.click();
    document.body.removeChild(elem);
  });

  document.getElementById("downloadBtn").addEventListener("click", () => {
    socket.emit("download:request", {});
  });

  document
    .getElementById("downloadSelectedBtn")
    .addEventListener("click", () => {
      socket.emit("download:request", {
        selection: world.getSelection(),
      });
    });
}
