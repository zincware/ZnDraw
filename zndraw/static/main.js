import { Cache } from "./pycom/Cache.js";
import { World } from "./World/World.js";
import { setUIEvents } from "./UI/UI.js";
import { initJSONEditor } from "./UI/json_editor.js";

function setupSocket() {
  const socket = io();
  socket.on("connect", function () {
    socket.emit("connection", { data: "I'm connected!" });
  });
  return socket;
}

function main() {
  const socket = setupSocket();
  const cache = new Cache(socket);
  const container = document.querySelector("#scene-container");
  const world = new World(container, cache, socket);

  setUIEvents(socket, world);
  initJSONEditor(socket);
  world.start();

  socket.emit("atoms:request", null, () => {
    world.setStep(0);
  });

  socket.on("display", (data) => {
    socket.emit("config");
    const slider = document.getElementById("frame-slider");
    slider.value = data.index;
    world.setStep(data.index);
  });

  // disable loading screen
  document.getElementById("atom-spinner").style.display = "none";
}

main();
