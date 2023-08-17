import { Cache } from "./pycom/Cache.js";



var socket = io();
socket.on("connect", function () {
    socket.emit("connection", { data: "I'm connected!" });
});

socket.emit("config", (data) => {
    const progress = document.getElementById("progress-visualizer");
    progress.ariaValueMax = data["n_frames"];
})

const cache = new Cache(socket);



function progressController(element) {
    document.getElementById(element.id).onpointerdown = function (event) {

        const controlls = document.getElementById("progress-controller");
        const progress = document.getElementById("progress-visualizer");
        const progressbarvisualizer = document.getElementById("progress-bar-visualizer");
        
        // TODO: set the hight once via css in the beginning
        controlls.style.top = event.clientY - controlls.offsetHeight / 2 + "px";
        
        function onPointerMove(event) {
            controlls.style.left = event.clientX - controlls.offsetWidth / 2 + "px";
            
            let percentage = (event.clientX / element.offsetWidth) * 100;
            progress.ariaValueNow = percentage / 100 * progress.ariaValueMax;
            progressbarvisualizer.style.width = percentage + "%";
            // progress.style.width = percentage + "%";
            
        }
        onPointerMove(event);

        document.onpointermove = onPointerMove;
        document.onpointerup = async function () {
            
            const atoms = await cache.get(Math.round(progress.ariaValueNow));
            console.log(atoms);
            

            document.onpointermove = null;
            document.onpointerup = null;
        }
    }
}

progressController(document.getElementById("progress-visualizer"));
