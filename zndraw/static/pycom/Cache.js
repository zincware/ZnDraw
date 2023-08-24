// Interface for the communication with Python to retrieve atoms

class Atom {
  constructor({ position, number, color, id, radius }) {
    this.position = position;
    this.number = number;
    this.id = id;
    this.color = color;
    this.radius = radius;
  }
}

class Atoms {
  constructor({ positions, cell, numbers, colors, radii, connectivity }) {
    this.positions = positions;
    this.cell = cell;
    this.numbers = numbers;
    this.colors = colors;
    this.radii = radii;
    this.connectivity = connectivity;

    this.length = this.positions.length;
  }

  [Symbol.iterator]() {
    let index = 0;
    const atoms = this; // this does not work
    return {
      next() {
        if (index < atoms.positions.length) {
          const atom = new Atom({
            position: atoms.positions[index],
            number: atoms.numbers[index],
            color: atoms.colors[index],
            radius: atoms.radii[index],
            id: index,
          });
          index++;
          return { value: atom, done: false };
        } else {
          return { done: true };
        }
      },
    };
  }
}

class Cache {
  constructor(socket) {
    this._socket = socket;
    this._cache = {};

    this._last_request = -999999;

    this._socket.on("cache:load", (data) => {
      this._get(data.index);
    });
    this._socket.on("cache:reset", () => {
      this.reset();
    });
  }

  async _get(id) {
    // if the id is not in the cache, request it from the server
    if (!(id in this._cache)) {
      // create a promise that resolves when the server response is received
      const dataPromise = new Promise((resolve) => {
        this._socket.emit("configuration:id", { id: id }, (data) => {
          Object.keys(data).forEach((key) => {
            this._cache[key] = new Atoms({
              positions: data[key].positions,
              cell: data[key].cell,
              numbers: data[key].numbers,
              colors: data[key].colors,
              radii: data[key].radii,
              connectivity: data[key].connectivity,
            });
          });
          resolve();
        });
      });
      // wait for the server response before returning the cached data
      await dataPromise;
    }
    return this._cache[id];
  }

  get(id) {
    // convert id to integer
    id = parseInt(id);
    const value = this._cache[id];

    if (value === undefined) {
      document.getElementById("frame-slider").disabled = true;

      if (Math.abs(id - this._last_request) < 15) {
        return value;
      }

      this._last_request = id;
      this._get(id);
    } else {
      document.getElementById("frame-slider").disabled = false;
      // set focus to the slider
      document.getElementById("frame-slider").focus();

      if (Math.abs(id - this._last_request) < 15) {
        return value;
      }

      for (let i = id + 1; i < id + 10; i++) {
        if (!(i in this._cache)) {
          this._last_request = i;
          this._get(i);
          break;
        }
      }
    }
    return value;
  }

  reset() {
    console.log("reset cache");
    this._cache = {};
    this._last_request = -999999;
  }
}

export { Cache };
