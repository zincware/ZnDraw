// Interface for the communication with Python to retrieve atoms
import * as THREE from "three";
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
  constructor({
    positions,
    cell,
    numbers,
    colors,
    radii,
    connectivity,
    calc,
    pbc,
  }) {
    this.positions = positions;
    this.cell = cell;
    this.numbers = numbers;
    this.colors = colors;
    this.radii = radii;
    this.connectivity = connectivity;
    this.calc = calc;
    this.pbc = pbc;

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
  select(indices) {
    const selectedPositions = indices.map((index) => this.positions[index]);
    const selectedNumbers = indices.map((index) => this.numbers[index]);
    const selectedColors = indices.map((index) => this.colors[index]);
    const selectedRadii = indices.map((index) => this.radii[index]);
    const selectedAtoms = new Atoms({
      positions: selectedPositions,
      cell: this.cell,
      numbers: selectedNumbers,
      colors: selectedColors,
      radii: selectedRadii,
      connectivity: this.connectivity,
      calc: this.calc,
      pbc: this.pbc,
    });
    return selectedAtoms;
  }

  getCenter() {
    const sum = this.positions.reduce((acc, position) => {
      const vec = new THREE.Vector3().fromArray(position);
      return acc.add(vec);
    }, new THREE.Vector3());

    const mean = sum.divideScalar(this.positions.length);

    return mean;
  }
}

class Cache {
  constructor(socket) {
    this._socket = socket;
    this._cache = {};
    this.world;

    this._socket.on("atoms:upload", (data) => {
      console.log(new Date().toISOString(), "Received atoms from Python");
      document.getElementById(
        "interaction-json-editor-submit",
      ).disabled = false;
      Object.keys(data).forEach((key) => {
        this._cache[key] = new Atoms({
          positions: data[key].positions,
          cell: data[key].cell,
          numbers: data[key].numbers,
          colors: data[key].colors,
          radii: data[key].radii,
          connectivity: data[key].connectivity,
          calc: data[key].calc,
          pbc: data[key].pbc,
        });
      });
      const slider = document.getElementById("frame-slider");
      slider.max = Object.keys(this._cache).length - 1;
      document.getElementById(
        "info",
      ).innerHTML = `${slider.value} / ${slider.max}`;
    });

    this._socket.on("atoms:delete", (ids) => {
      for (const id of ids) {
        delete this._cache[id];
      }
      // move all keys after id one step back
      const remainingKeys = Object.keys(this._cache);
      for (let i = ids[0]; i < remainingKeys.length; i++) {
        const currentKey = remainingKeys[i];
        const newIndex = i;
        if (currentKey !== newIndex) {
          this._cache[newIndex] = this._cache[currentKey];
          delete this._cache[currentKey];
        }
      }
      // update slider
      const slider = document.getElementById("frame-slider");

      // update world
      if (this.world.getStep() >= slider.max) {
        this.world.setStep(remainingKeys.length - 1);
      } else {
        let newStep = this.world.getStep();
        ids.forEach((id) => {
          if (this.world.getStep() > id) {
            newStep = newStep - 1;
          }
        });
        this.world.setStep(newStep);
      }

      slider.max = remainingKeys.length - 1;
      document.getElementById(
        "info",
      ).innerHTML = `${slider.value} / ${slider.max}`;
    });

    this._socket.on("atoms:insert", (data) => {
      // move all keys after id one step forward
      const remainingKeys = Object.keys(this._cache);
      const id = parseInt(Object.keys(data)[0]);
      for (let i = remainingKeys.length - 1; i >= id; i--) {
        const currentKey = remainingKeys[i];
        const newIndex = i + 1;
        this._cache[newIndex] = this._cache[currentKey];
        delete this._cache[currentKey];
      }
      // insert new atoms
      this._cache[id] = new Atoms({
        positions: data[id].positions,
        cell: data[id].cell,
        numbers: data[id].numbers,
        colors: data[id].colors,
        radii: data[id].radii,
        connectivity: data[id].connectivity,
        calc: data[id].calc,
        pbc: data[id].pbc,
      });
      // update slider
      const slider = document.getElementById("frame-slider");
      slider.max = Object.keys(this._cache).length - 1;
      document.getElementById(
        "info",
      ).innerHTML = `${slider.value} / ${slider.max}`;
    });

    this._socket.on(
      "atoms:download",
      function (ids, callback) {
        // send all atoms at once
        const data = {};
        ids.forEach((x) => {
          data[x] = this._cache[x];
        });
        callback(data);
      }.bind(this),
    );

    this._socket.on("atoms:clear", (start_index) => {
      // remove everything from the cache starting from start_index
      Object.keys(this._cache).forEach((key) => {
        if (parseInt(key) >= start_index) {
          delete this._cache[key];
        }
      });
    });

    this._socket.on(
      "atoms:length",
      function (callback) {
        callback(this.get_length());
      }.bind(this),
    );
  }

  get(id) {
    // convert id to integer
    id = parseInt(id);
    return this._cache[id];
  }

  get_length() {
    return Object.keys(this._cache).length;
  }

  getAllAtoms() {
    return this._cache;
  }

  attachWorld(world) {
    this.world = world;
  }
}

export { Cache };
