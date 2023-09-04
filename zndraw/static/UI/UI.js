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

function setupNavbarLeft() {

  function showMenu(menu) {
    const menus = ["selectionMenu", "interactionMenu", "sceneMenu", "drawMenu", "analysisMenu"];
    for (let i = 0; i < menus.length; i++) {
      if ((menus[i] === menu) && (document.getElementById(menus[i]).style.display === "none")) {
        document.getElementById(menus[i]).style.display = "block";
        document.getElementById(`${menus[i]}Btn`).classList.add("active");
      } else {
        document.getElementById(menus[i]).style.display = "none";
        document.getElementById(`${menus[i]}Btn`).classList.remove("active");
      }
    }
  }

  function closeMenu(menu) {
    document.getElementById(menu).style.display = "none";
    document.getElementById(`${menu}Btn`).classList.remove("active");
  }

  document.getElementById("selectionMenuBtn").onclick = () => {
    showMenu("selectionMenu");
  };
  document.getElementById("selectionMenuClose").onclick = () => {
    closeMenu("selectionMenu");
  };

  document.getElementById("interactionMenuBtn").onclick = () => {
    showMenu("interactionMenu");
  };
  document.getElementById("interactionMenuClose").onclick = () => {
    closeMenu("interactionMenu");
  };

  document.getElementById("sceneMenuBtn").onclick = () => {
    showMenu("sceneMenu");
  };
  document.getElementById("sceneMenuClose").onclick = () => {
    closeMenu("sceneMenu");
  };

  document.getElementById("drawMenuBtn").onclick = () => {
    showMenu("drawMenu");
  };
  document.getElementById("drawMenuClose").onclick = () => {
    closeMenu("drawMenu");
  };

  document.getElementById("analysisMenuBtn").onclick = () => {
    showMenu("analysisMenu");
  };
  document.getElementById("analysisMenuClose").onclick = () => {
    closeMenu("analysisMenu");
  };
}

export function setUIEvents(socket, cache, world) {
  resizeOffcanvas();
  setupUpload(socket);
  setupNavbarLeft();

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
