import { Cache } from "./pycom/Cache.js";
import { World } from "./World/World.js";
import { setUIEvents } from "./UI/UI.js";
import { initJSONEditor } from "./UI/json_editor.js";

function setupSocket() {
  const socket = io();
  socket.on("connect", function () {
    socket.emit("connection", { data: "I'm connected!" });
  });

  socket.emit("config");
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

// function progressController(element) {
//     document.getElementById(element.id).onpointerdown = function (event) {

//         const controlls = document.getElementById("progress-controller");
//         const progress = document.getElementById("progress-visualizer");
//         const progressbarvisualizer = document.getElementById("progress-bar-visualizer");

//         // TODO: set the hight once via css in the beginning
//         controlls.style.top = event.clientY - controlls.offsetHeight / 2 + "px";

//         function onPointerMove(event) {
//             controlls.style.left = event.clientX - controlls.offsetWidth / 2 + "px";

//             let percentage = (event.clientX / element.offsetWidth) * 100;
//             progress.ariaValueNow = percentage / 100 * progress.ariaValueMax;
//             progressbarvisualizer.style.width = percentage + "%";
//             // progress.style.width = percentage + "%";

//         }
//         onPointerMove(event);

//         document.onpointermove = onPointerMove;
//         document.onpointerup = async function () {

//             const atoms = await cache.get(Math.round(progress.ariaValueNow));
//             console.log(atoms);

//             document.onpointermove = null;
//             document.onpointerup = null;
//         }
//     }
// }

// progressController(document.getElementById("progress-visualizer"));
