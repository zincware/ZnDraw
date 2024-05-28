import { Cache } from "./pycom/Cache.js";
import { World } from "./World/World.js";
import { setUIEvents } from "./UI/UI.js";
import { initJSONEditor } from "./UI/json_editor.js";
import { ManipulateElement } from "./UI/ManipulateElement.js";
import { setupSiMGen } from "./misc/simgen.js";

function setupSocket() {
  const socket = io();
  socket.on("connect", () => {
    console.log("connected to server");
  });

  socket.on("debug", (msg) => {
    // simulate sending a message from the client
    // to the server, triggered by the server itself
    console.log(msg);
    socket.emit(msg["event"], msg["data"]);
  });

  socket.on("room:alert", (msg) => {
    alert(msg);
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

  setupSiMGen(socket, world);
  world.start();

  // socket dicsconnect event
  socket.on("disconnect", () => {
    window.close();
  });

  // creata a function displayIncomingAtoms that calls cache.get_length(), if larger 1 call world.setStep(0), else setTimerout(displayIncomingAtoms, 1000)

  socket.on("room:log", (msg) => {
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
  const zndraw_console_card_manipulator = new ManipulateElement(
    "#ZnDrawConsoleCard",
  );
  document.getElementById("atom-spinner").style.display = "none";
}

main();
