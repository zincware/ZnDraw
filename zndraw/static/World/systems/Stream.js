import { Clock, SkeletonHelper } from 'three';

const clock = new Clock();

// const socket = 

// // Connection opened
// socket.addEventListener("open", (event) => {
//   socket.send("Hello Server!");
// });

let socket;


class Stream {
  constructor() {
    this._ready = false;
    this.data = null;
    this.step = 0;
    socket = new WebSocket('ws://' + location.host + '/echo');

    socket.onopen = (event) => {
      console.log("Websocket connection opened");
      this._ready = true;
    };
    socket.addEventListener("message", (event) => {
      this._ready = true;
      this.data = JSON.parse(event.data);
    });
  }

  requestFrame() {
    if (this._ready) {
      this._ready = false;
      socket.send(JSON.stringify({ "step": this.step }));
    }
  }

  get_next_frame() {
    console.log("Step" + this.step);
    if (this.data == null) {
      return undefined;
    }
    try {
      let data = this.data[this.step];
      delete this.data[this.step];
      if (data !== undefined) {
        this.step += 1;
      } else {
        this.step = 0;
      }
      return data;
    } catch (error) {
      return undefined;
    }
  }
}

export { Stream };
