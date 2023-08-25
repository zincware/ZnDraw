import { Cache } from "./pycom/Cache.js";
import { World, Player } from "./World/World.js";
import { setUIEvents } from "./UI/UI.js";
import { initJSONEditor } from "./UI/json_editor.js";

function setupSocket() {
  const socket = io();
  socket.on("connect", function () {
    socket.emit("connection", { data: "I'm connected!" });
  });
  return socket;
}


function setupSlider(socket, world) {
  
  
}

function main() {
  const socket = setupSocket();
  const cache = new Cache(socket);
  const container = document.querySelector("#scene-container");
  const world = new World(container, cache, socket);
  const player = new Player(world, cache);

  setUIEvents(socket, world);
  initJSONEditor(socket, cache, world);
  world.start();

  socket.emit("atoms:request", null, () => {
    world.setStep(0);    
  });

  

  // disable loading screen
  document.getElementById("atom-spinner").style.display = "none";
}

main();
