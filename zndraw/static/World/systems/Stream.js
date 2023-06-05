class Stream {
  constructor(config) {
    this.step = 1;
    this.config = config;
    this.data = null;
    this.last_request = 1;
    this._buffer_filled = false;

    // fetch load to start loading data in the background
    fetch("/load");

    this.setup_event_source();
  }

  setup_event_source() {
    this.eventSource = new EventSource("/frame-stream");

    this.eventSource.onmessage = (event) => {
      const data = JSON.parse(event.data);
      // if data length is zero, close the connection
      if (Object.keys(data).length === 0) {
        this.eventSource.close();
        return;
      }
      this.data = { ...this.data, ...data };
      for (const key in this.data) {
        if (key < this.step - this.config.config.js_frame_buffer[0]) {
          delete this.data[key];
        }
        if (key > this.step + this.config.config.js_frame_buffer[1]) {
          delete this.data[key];
        }
      }
    };
  }

  requestFrame() {
    // fetch frame-set with post request step: this.step
    if (this.step === this.last_request) {
      return;
    }
    this.last_request = this.step;
    console.log("Requesting frame " + this.step);
    fetch("/frame-set", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({ step: this.step }),
    }).then((response) => {
      // close the event source and reopen it to get the new data
      this.eventSource.close();
      this.setup_event_source();
    });
  }

  deleteCache() {
    this.data = null;
    this.last_request = -1;
    this.requestFrame();
  }

  setStep() {
    this.step = step;
    this.requestFrame();
  }

  getStep() {
    return this.step;
  }

  get_next_frame() {
    if (this.data == null) {
      return undefined;
    }
    console.log(
      "Step " +
        this.step +
        " with cache size: " +
        Object.keys(this.data).length,
    );
    const data = this.data[this.step];
    if (data !== undefined) {
      // TODO this also happens if the stream is to slow to keep up!
      if (this.config.play) {
        this.set_step(this.step + 1);
      }
    }
    // TODO: do we need to request new fames, if they are already in the buffer?
    if (
      this.step - this.last_request >
        this.config.config.js_frame_buffer[1] / 2 ||
      this.last_request - this.step > this.config.config.js_frame_buffer[0] / 2
    ) {
      this.requestFrame();
    }
    if (
      this.step > this.config.config.total_frames &&
      this.config.config.total_frames > 0
    ) {
      // temporary freeze for larger than 1000
      // TODO we need a good way for handling jumps in frames
      if (this.config.config.auto_restart) {
        this.set_step(0);
      } else {
        this.set_step(Math.max(0, this.step - 1));
      }
    }
    return data;
  }
}

export { Stream };
