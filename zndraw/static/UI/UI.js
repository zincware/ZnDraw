function resizeOffcanvas() {
  // Rescale offcanvas by dragging
  let active_offcanvas_border;
  const offcanvas_borders = document.getElementsByClassName("offcanvas-border");

  function resize_offcanvas(e) {
    if (e.clientX < 200) {
      return;
    }
    active_offcanvas_border.parentNode.style.width = `${e.clientX}px`;
    active_offcanvas_border.style.left = `${e.clientX}px`;
  }

  for (let i = 0; i < offcanvas_borders.length; i++) {
    offcanvas_borders[i].onpointerdown = function (e) {
      active_offcanvas_border = this;
      document.addEventListener("pointermove", resize_offcanvas);
    };
  }

  document.addEventListener("pointerup", (e) => {
    document.removeEventListener("pointermove", resize_offcanvas);
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

export function setUIEvents(socket, cache, world) {
  resizeOffcanvas();
  setupUpload(socket);

  document.getElementById("ExitBtn").addEventListener("click", () => {
    fetch("/exit", { method: "GET" });
  });

  document.getElementById("downloadBtn").addEventListener("click", () => {
    socket.emit("download", { atoms_list: cache.getAllAtoms() }, (data) => {
      const blob = new Blob([data], { type: "text/csv" });
      const elem = window.document.createElement("a");
      elem.href = window.URL.createObjectURL(blob);
      elem.download = "trajectory.xyz";
      document.body.appendChild(elem);
      elem.click();
      document.body.removeChild(elem);
    });
  });

  document
    .getElementById("downloadSelectedBtn")
    .addEventListener("click", () => {
      socket.emit(
        "download",
        {
          atoms_list: cache.getAllAtoms(),
          selection: world.getSelection(),
        },
        (data) => {
          const blob = new Blob([data], { type: "text/csv" });
          const elem = window.document.createElement("a");
          elem.href = window.URL.createObjectURL(blob);
          elem.download = "trajectory.xyz";
          document.body.appendChild(elem);
          elem.click();
          document.body.removeChild(elem);
        },
      );
    });

  const helpBtn = document.getElementById("HelpBtn");

  // helpBtn.addEventListener("mouseover", () => {
  //   new bootstrap.Collapse(document.getElementById("helpBoxCollapse"), {
  //     toggle: false,
  //   }).show();
  // });
  // helpBtn.addEventListener("mouseout", () => {
  //   new bootstrap.Collapse(document.getElementById("helpBoxCollapse"), {
  //     toggle: false,
  //   }).hide();
  // });
}
