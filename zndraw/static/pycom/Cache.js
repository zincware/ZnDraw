// Interface for the communication with Python to retrieve atoms

class Atom {
  constructor({
    position, number, color, id, radius,
  }) {
    this.position = position;
    this.number = number;
    this.id = id;
    this.color = color;
    this.radius = radius;
  }
}

class Atoms {
  constructor({
    positions, cell, numbers, colors, radii, connectivity, calc,
  }) {
    this.positions = positions;
    this.cell = cell;
    this.numbers = numbers;
    this.colors = colors;
    this.radii = radii;
    this.connectivity = connectivity;
    this.calc = calc;

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
        }
        return { done: true };
      },
    };
  }
}

class Cache {
  constructor(socket) {
    this._socket = socket;
    this._cache = {};

    this._last_request = -999999;

    this._socket.on('atoms:upload', (data) => {
      console.log('Received atoms from Python');
      document.getElementById('interaction-json-editor-submit').disabled = false;
      Object.keys(data).forEach((key) => {
        this._cache[key] = new Atoms({
          positions: data[key].positions,
          cell: data[key].cell,
          numbers: data[key].numbers,
          colors: data[key].colors,
          radii: data[key].radii,
          connectivity: data[key].connectivity,
          calc: data[key].calc,
        });
      });
      const slider = document.getElementById('frame-slider');
      slider.max = Object.keys(this._cache).length - 1;
      document.getElementById(
        'info',
      ).innerHTML = `${slider.value} / ${slider.max}`;
    });

    this._socket.on('atoms:download', (ids) => {
      // iterate through the ids and emit the atoms object for each
      // ids.forEach((x) => {
      //   // let data = {};
      //   // data[x] = this._cache[x];
      //   // this._socket.emit('atoms:download', data);
      //   this._socket.emit('atoms:download', this._cache[x]);
      // });

      // send all atoms at once
      const data = {};
      ids.forEach((x) => {
        data[x] = this._cache[x];
      });
      this._socket.emit('atoms:download', data);
    });

    this._socket.on('atoms:clear', (start_index) => {
      // remove everything from the cache starting from start_index
      Object.keys(this._cache).forEach((key) => {
        if (parseInt(key) >= start_index) {
          delete this._cache[key];
        }
      });
    });

    this._socket.on('atoms:size', () => {
      this._socket.emit('atoms:size', this.get_length());
    });
  }

  get(id) {
    // convert id to integer
    id = parseInt(id);
    return this._cache[id];
  }

  reset() {
    console.log('reset cache');
    this._cache = {};
    this._last_request = -999999;
  }

  get_length() {
    return Object.keys(this._cache).length;
  }

  getAllAtoms() {
    return this._cache;
  }
}

export { Cache };
