import { Cache } from "./pycom/Cache.js";
import { World } from "./World/World.js";
import { setUIEvents } from "./UI/UI.js";
import { initJSONEditor } from "./UI/json_editor.js";

function setupSocket() {
  const socket = io();
  socket.on("connect", () => {
    console.log("connected to server");
  });

  socket.on("connect_error", (err) => {
    console.log("connection could not be established - trying again.");
    // try to connect again
    setTimeout(() => {
      if (socket.connected) {
        return;
      }
      socket.connect();
    }, 100);
  });
  // if socket.on
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

  socket.on("message:log", (msg) => {
    console.log(msg);
    // add <p> to console body
    const consoleBody = document.getElementById("ZnDrawConsoleBody");
    const p = document.createElement("p");
    // ISO 8601 format + msg
    p.innerHTML = new Date().toISOString() + " " + msg;
    consoleBody.appendChild(p);
    // automatically scroll to bottom of consoleBody
    consoleBody.scrollTop = consoleBody.scrollHeight;
  });

  document.getElementById("atom-spinner").style.display = "none";
}

main();
