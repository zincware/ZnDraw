import { World } from "./World/World.js";
import { setUIEvents } from "./UI/UI.js";

// create a config class that repeatedly calls the config endpoint and updates the config
class Config {
  constructor() {
    this.config = {};
    this.config_url = "/config";
  }

  start() {
    this.update_config();
  }

  update_config() {
    fetch(this.config_url)
      .then((response) => response.json())
      .then((data) => {
        this.config = data;
        if (this.onLoadCallback) {
          this.onLoadCallback(this);
          this.onLoadCallback = null;
        }
        setTimeout(() => this.update_config(), 100);
      });
  }

  onLoad(callback) {
    this.onLoadCallback = callback;
  }

  update(config) {
    // this.config = {...this.config, ...config};
    fetch(this.config_url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(config),
    });
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
