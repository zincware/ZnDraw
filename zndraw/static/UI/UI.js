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
    const menus = [
      "selectionMenu",
      "interactionMenu",
      "sceneMenu",
      "drawMenu",
      "analysisMenu",
    ];
    for (let i = 0; i < menus.length; i++) {
      if (
        menus[i] === menu &&
        document.getElementById(menus[i]).style.display === "none"
      ) {
        document.getElementById(menus[i]).style.display = "block";
        document.getElementById(`${menus[i]}Btn`).classList.add("active");
      } else {
        document.getElementById(menus[i]).style.display = "none";
        document.getElementById(`${menus[i]}Btn`).classList.remove("active");
      }
    }
  }

  const popovers = {
    selectionMenu: new bootstrap.Popover(
      document.getElementById("selectionMenuBtn"),
    ),
    interactionMenu: new bootstrap.Popover(
      document.getElementById("interactionMenuBtn"),
    ),
    sceneMenu: new bootstrap.Popover(document.getElementById("sceneMenuBtn")),
    drawMenu: new bootstrap.Popover(document.getElementById("drawMenuBtn")),
    analysisMenu: new bootstrap.Popover(
      document.getElementById("analysisMenuBtn"),
    ),
  };

  function closeMenu(menu) {
    document.getElementById(menu).style.display = "none";
    document.getElementById(`${menu}Btn`).classList.remove("active");
  }

  // close all popovers when clicking anywhere
  document.addEventListener("click", () => {
    popovers.selectionMenu.hide();
    popovers.interactionMenu.hide();
    popovers.sceneMenu.hide();
    popovers.drawMenu.hide();
    popovers.analysisMenu.hide();
  });

  document.getElementById("selectionMenuBtn").onclick = () => {
    showMenu("selectionMenu");
    popovers.selectionMenu.hide();
  };

  document.getElementById("selectionMenuClose").onclick = () => {
    closeMenu("selectionMenu");
  };

  document.getElementById("interactionMenuBtn").onclick = () => {
    showMenu("interactionMenu");
    popovers.interactionMenu.hide();
  };
  document.getElementById("interactionMenuClose").onclick = () => {
    closeMenu("interactionMenu");
  };

  document.getElementById("sceneMenuBtn").onclick = () => {
    showMenu("sceneMenu");
    popovers.sceneMenu.hide();
  };
  document.getElementById("sceneMenuClose").onclick = () => {
    closeMenu("sceneMenu");
  };

  document.getElementById("drawMenuBtn").onclick = () => {
    showMenu("drawMenu");
    popovers.drawMenu.hide();
  };
  document.getElementById("drawMenuClose").onclick = () => {
    closeMenu("drawMenu");
  };

  document.getElementById("analysisMenuBtn").onclick = () => {
    showMenu("analysisMenu");
    popovers.analysisMenu.hide();
  };
  document.getElementById("analysisMenuClose").onclick = () => {
    closeMenu("analysisMenu");
  };

  document.getElementById("drawMenuBtn").onpointerenter = () => {
    if (document.getElementById("drawMenu").style.display === "none") {
      popovers.drawMenu.show();
    }
  };
  document.getElementById("drawMenuBtn").onpointerleave = () => {
    popovers.drawMenu.hide();
  };

  document.getElementById("sceneMenuBtn").onpointerenter = () => {
    if (document.getElementById("sceneMenu").style.display === "none") {
      popovers.sceneMenu.show();
    }
  };
  document.getElementById("sceneMenuBtn").onpointerleave = () => {
    popovers.sceneMenu.hide();
  };

  document.getElementById("selectionMenuBtn").onpointerenter = () => {
    if (document.getElementById("selectionMenu").style.display === "none") {
      popovers.selectionMenu.show();
    }
  };
  document.getElementById("selectionMenuBtn").onpointerleave = () => {
    popovers.selectionMenu.hide();
  };

  document.getElementById("interactionMenuBtn").onpointerenter = () => {
    if (document.getElementById("interactionMenu").style.display === "none") {
      popovers.interactionMenu.show();
    }
  };
  document.getElementById("interactionMenuBtn").onpointerleave = () => {
    popovers.interactionMenu.hide();
  };

  document.getElementById("analysisMenuBtn").onpointerenter = () => {
    if (document.getElementById("analysisMenu").style.display === "none") {
      popovers.analysisMenu.show();
    }
  };
  document.getElementById("analysisMenuBtn").onpointerleave = () => {
    popovers.analysisMenu.hide();
  };
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

export function setUIEvents(socket, cache, world) {
  // resizeOffcanvas();
  setupUpload(socket);
  setupNavbarLeft();
  setupMobile();

  document.getElementById("ExitBtn").addEventListener("click", () => {
    fetch("/exit", { method: "GET" });
  });

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
