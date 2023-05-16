import { World } from "./World/World.js";
import { setUIEvents } from "./UI/UI.js";

// create a config class that repeatedly calls the config endpoint and updates the config
class Config {
  constructor() {
    this.config = {};
    this.config_url = "/config";
    this.step = 0;
    this.set_step_callbacks = [];
    this.rebuild_callbacks = [];
    this.play = true;
    this.selected = [];
    this.pressed_keys = {};
  }

  start() {
    this.update_config();
  }

  rebuild() {
    console.log("rebuild");
    for (const callback of this.rebuild_callbacks) {
      this.update_config(0).then(() => callback(this));
    }
  }

  async update_config(timeout = 100) {
    await fetch(this.config_url)
      .then((response) => response.json())
      .then((data) => {
        this.config = data;
        if (this.onLoadCallback) {
          this.onLoadCallback(this);
          this.onLoadCallback = null;
        }
        if (timeout > 0) {
          setTimeout(() => this.update_config(), timeout);
        }
      });
  }

  onLoad(callback) {
    this.onLoadCallback = callback;
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
  }
}

function main() {
  // Get a reference to the container element
  const container = document.querySelector("#scene-container");
  const config = new Config();
  config.onLoad(setUIEvents);
  config.start();

  // 1. Create an instance of the World app
  const world = new World(container, config);

  // 2. Render the scene
  world.start();
}

main();
