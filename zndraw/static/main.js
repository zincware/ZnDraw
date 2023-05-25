import { World } from "./World/World.js";
import { setUIEvents } from "./UI/UI.js";

// create a config class that repeatedly calls the config endpoint and updates the config
class Config {
  constructor() {
    this.config = {};
    this.config_url = "/config";
    this.step = 0;
    this.set_step_callbacks = [];
    this.play = true;
    this.selected = [];
    this.pressed_keys = {};
    this.draw_vectors = [];
  }

  start() {
    this.update_config();
  }

  async update_config(timeout = 100) {
    await fetch(this.config_url)
      .then((response) => response.json())
      .then((data) => {
        this.config = data;
        if (this.onLoadCallback) {
          this.onLoadCallback();
          this.onLoadCallback = null;
        }
        document.getElementById(
          "info",
        ).innerHTML = `${this.step} / ${this.config.total_frames}`;
        if (timeout > 0) {
          setTimeout(() => this.update_config(), timeout);
        }
      });
  }

  async update(config) {
    // this.config = {...this.config, ...config};
    await fetch(this.config_url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(config),
    });
  }

  set_step(step) {
    this.step = step;
    for (const callback of this.set_step_callbacks) {
      callback(this);
    }
    document.getElementById(
      "info",
    ).innerHTML = `${this.step} / ${this.config.total_frames}`;
  }
}

function main() {
  // Get a reference to the container element
  const container = document.querySelector("#scene-container");
  const config = new Config();

  // 1. Create an instance of the World app
  const world = new World(container, config);

  // config.onLoad(setUIEvents);
  config.onLoadCallback = () => {
    setUIEvents(config, world);
  };
  config.start();

  // 2. Render the scene
  world.start();
}

main();
