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

    this._socket.on("atoms:upload", (data) => {
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

      const slider = document.getElementById("frame-slider");
      slider.max = Object.keys(this._cache).length - 1;
    });

    this._socket.on("atoms:download", (ids) => {
      // iterate through the ids and emit the atoms object for each
      ids.forEach((id) => {
        this._socket.emit("atoms:download", this._cache[id]);
      });
    });

    this._socket.on("atoms:clear", (start_index) => {
      // remove everything from the cache starting from start_index
      Object.keys(this._cache).forEach((key) => {
        if (parseInt(key) >= start_index) {
          delete this._cache[key];
        }
      });
    });
  }

  get(id) {
    // convert id to integer
    id = parseInt(id);
    return this._cache[id];
  }

  reset() {
    console.log("reset cache");
    this._cache = {};
    this._last_request = -999999;
  }

  get_length() {
    return Object.keys(this._cache).length;
  }
}

export { Cache };
