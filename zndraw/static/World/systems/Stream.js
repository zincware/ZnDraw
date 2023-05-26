class Stream {
  constructor(config) {
    this.config = config;
    this.step = 1;
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
        if (key < this.step - 10) {
          delete this.data[key];
        }
        if (key > this.step + 100) {
          delete this.data[key];
        }
      }
    };
  }

  requestFrame() {
    // fetch frame-set with post request step: this.config.step
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
      // if the event source is closed, open it again
      if (this.eventSource.readyState === 2) {
        this.setup_event_source();
      }
    });
  }

  deleteCache() {
    this.data = null;
    this.last_request = -1;
    this.requestFrame();
  }

  setStep(step) {
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
        this.setStep(this.step + 1);
      }
    } else {
      // if the data is not available, request it. This should not happen if the stream is fast enough
      console.log("Unexpected request frame " + this.step);
      this.requestFrame();
    }
    if (
      this.step - this.last_request > 50 ||
      this.step < this.last_request // Going backwards is stil a problem
    ) {
      this.requestFrame();
    }
    if (this.step >= this.config.total_frames) {
      // temporary freeze for larger than 1000
      // TODO we need a good way for handling jumps in frames
      this.config.step = 0;
    }
    return data;
  }
}

export { Stream };
