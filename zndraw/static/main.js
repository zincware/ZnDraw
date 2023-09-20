import { Cache } from "./pycom/Cache.js";
import { World } from "./World/World.js";
import { setUIEvents } from "./UI/UI.js";
import { initJSONEditor } from "./UI/json_editor.js";

function setupSocket() {
  const socket = io();
  socket.on("connect", () => {
    socket.emit("connection", { data: "I'm connected!" });
  });
  return socket;
}

function main() {
  const socket = setupSocket();
  const cache = new Cache(socket);
  const container = document.querySelector("#scene-container");
  const world = new World(container, cache, socket);

  setUIEvents(socket, cache, world);
  initJSONEditor(socket, cache, world);
  world.start();

  // socket dicsconnect event
  socket.on("disconnect", () => {
    window.close();
  });

  // creata a function displayIncomingAtoms that calls cache.get_length(), if larger 1 call world.setStep(0), else setTimerout(displayIncomingAtoms, 1000)

  const displayIncomingAtoms = () => {
    cache.get_length();
    if (cache.get_length() > 0) {
      world.setStep(0);
      document.getElementById("atom-spinner").style.display = "none";
    } else {
      setTimeout(displayIncomingAtoms, 100);
    }
  };

  socket.emit("atoms:request", window.location.href, () => {
    displayIncomingAtoms();
  });

  socket.on("message:log", (msg) => {
    console.log(msg);
  });

  // keypress J
  document.addEventListener("keydown", (event) => {
    if (event.key === "j") {
      socket.emit("numpy", (data) => {
        // create FLOAT32ARRAY from data
        const array = new Float64Array(data);
        console.log(array);     
      }
      );
    }
  });
}

main();
