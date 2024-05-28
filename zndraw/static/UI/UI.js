function setupInfo() {
  fetch("static/info.md")
    .then((response) => response.text())
    .then((text) => {
      document.getElementById("helpModalBody").innerHTML = marked.parse(text);
    });

  const tutorialIframe = document.getElementById("tutorialIframe");
  if (tutorialIframe) {
    // set width and height of iframe to 100% of parent
    tutorialIframe.style.width = "100%";
    // set height to 80% of screen height
    tutorialIframe.style.height = window.innerHeight * 0.7 + "px";
  }
}

function setupCodeAccess() {
  // in data-token is the token
  const token = document.getElementById("token").dataset.token;
  const url = window.location.href.replace(/\/$/, "");

  fetch("static/code.md")
    .then((response) => response.text())
    .then((text) => {
      // replace {{token}} with token
      text = text.replace("{{token}}", token);
      text = text.replace("{{url}}", url);
      document.getElementById("codeModalContent").innerHTML =
        marked.parse(text);
    });
}

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
    socket.emit("file:upload", {
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
    // Disabled, because not reliable in combination with three.js
    // mouse movements.

    // show the overlay as long as the file is dragged over the scene
    // const overlay = document.getElementById("overlay");
    // overlay.style.display = "block";
  });

  scene.addEventListener("dragleave", (event) => {
    event.preventDefault();
    // hide the overlay when the file is dragged out of the scene
    // const overlay = document.getElementById("overlay");
    // overlay.style.display = "none";
  });

  scene.addEventListener("drop", (event) => {
    event.preventDefault();
    // // hide the overlay when the file is dropped
    // const overlay = document.getElementById("overlay");
    // overlay.style.display = "none";

    // read the file
    const file = event.dataTransfer.files[0];
    if (!file) {
      console.error("No file was dropped");
      return;
    }

    const reader = new FileReader();
    reader.readAsBinaryString(file);

    // send the file to the server
    reader.addEventListener("load", () => {
      socket.emit("file:upload", {
        content: reader.result,
        filename: file.name,
      });
    });
  });
}

function setupTrashClick(socket) {
  document.getElementById("trashBtn").addEventListener("click", () => {
    socket.emit("modifier:run", {
      method: { discriminator: "ClearTools" },
    });
    setTimeout(function () {
      document.getElementById("drawRemoveCanvas").click();
    }, 1000);
    // TODO: should this be a modifier?
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

function setupFrameInput(world) {
  const frame_input = document.getElementById("frameInput");
  const frame_input_modal = document.getElementById("frameInputModal");

  frame_input_modal.addEventListener("shown.bs.modal", () => {
    frame_input.focus();
  });

  frame_input.addEventListener("change", () => {
    const modal = bootstrap.Modal.getInstance(frame_input_modal);
    modal.hide();

    if (
      parseInt(frame_input.value) >
        document.getElementById("frameProgress").ariaValueMax ||
      parseInt(frame_input.value) < 0
    ) {
      alert("The frame you entered is out of range.");
    } else {
      const step = parseInt(frame_input.value);
      world.setStep(step);
    }

    frame_input.value = "";
  });
}

function setupConnectedUsers(socket) {
  const dropdown = document.getElementById("connectedUsersDropdown");
  if (dropdown === null) {
    return;
  }
  // for each user connected user there is
  // row > col py-1 d-grid > btn btn-outline-secondary
  // row > col py-2 d-grid > form-check-input (step)
  // row > col py-2 d-grid > form-check-label (camera)

  const token = document.getElementById("token").dataset.token;
  const url = window.location.href.replace(/\/$/, "");
  const toastLiveExample = document.getElementById("liveToast");
  const toast = new bootstrap.Toast(toastLiveExample);
  const toastBody = document.getElementById("toastBody");

  const copyTokenUrlBtns = document.getElementsByClassName("copyTokenUrlBtn");
  Array.from(copyTokenUrlBtns).forEach((btn) => {
    btn.addEventListener("click", () => {
      navigator.clipboard.writeText(url + "/token/" + token);
      // show text including the URL that was copied
      toastBody.innerHTML =
        "Copied URL to clipboard: " + url + "/token/" + token;
      toast.show();
    });
  });

  let name;

  socket.on("room:users:refresh", (data) => {
    // remove all rows except the first one (header)
    const rows = dropdown.querySelectorAll(".row");

    // save the shared camera state
    const cameraConnectedUser = document.querySelector(
      'input[name="sharedCamera"]:checked',
    );

    const entries = [];
    rows.forEach((row, index) => {
      if (index > 0) {
        row.remove();
      }
    });
    data.forEach((username, idx) => {
      if (name === undefined) {
        name = username;
      }
      const row = document.createElement("div");
      row.classList.add("row");
      const col1 = document.createElement("div");
      col1.classList.add("col-6");
      col1.classList.add("py-1");
      col1.classList.add("d-grid");
      const btn = document.createElement("button");
      btn.classList.add("btn");
      btn.classList.add("btn-outline-secondary");
      btn.disabled = true;
      btn.innerHTML = username;
      col1.appendChild(btn);
      row.appendChild(col1);

      const col3 = document.createElement("div");
      col3.classList.add("col");
      col3.classList.add("py-2");
      const inputCamera = document.createElement("input");
      inputCamera.classList.add("form-check-input");
      inputCamera.type = "radio";
      inputCamera.name = "sharedCamera";
      inputCamera.id = inputCamera.name + "-" + idx;
      inputCamera.value = username;

      inputCamera.addEventListener("change", () => {
        socket.emit("connectedUsers:subscribe:camera", {
          user: username,
        });
      });

      col3.appendChild(inputCamera);
      row.appendChild(col3);

      if (username === name) {
        inputCamera.checked = true;
        btn.innerHTML = username + " (you)";
      }
      if (username === cameraConnectedUser?.value) {
        inputCamera.checked = true;
      }

      entries.push({ name: username, row: row });
    });
    // sort the entries such that the first is name==name
    entries.sort((a, b) => {
      if (a.name === name) {
        return -1;
      }
      if (b.name === name) {
        return 1;
      }
      return a.name.localeCompare(b.name);
    });
    entries.forEach((entry) => {
      dropdown.appendChild(entry.row);
    });
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
  setupFrameInput(world);
  setupConnectedUsers(socket);
  setupInfo();
  setupCodeAccess();

  socket.on("file:download", (data) => {
    const blob = new Blob([data], { type: "text/csv" });
    const elem = window.document.createElement("a");
    elem.href = window.URL.createObjectURL(blob);
    elem.download = "trajectory.xyz";
    document.body.appendChild(elem);
    elem.click();
    document.body.removeChild(elem);
  });

  document.getElementById("downloadBtn").addEventListener("click", () => {
    socket.emit("file:download");
  });
}
