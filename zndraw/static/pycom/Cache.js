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
    this._cache = {};
    this.world;
    this.socket = socket;
    this._requested = {};
    this._size_timeout = false;
    this._received_data = false;

    this._length = 0;
    this.socket.on("room:frames:refresh", (ids) => {
      ids = ids.map((id) => parseInt(id));
      console.log("Cache refresh", ids);
      this.socket.emit("room:frames:get", ids, (data) => {
        this.setFrames(data);
        const currentStep = this.world.getStep();
        if (ids.includes(currentStep)) {
          this.world.setStep(currentStep);
        }
      });
    });
    // call get_length periodically
    setInterval(() => {
      this.get_length();
    }, 1000);
  }

  get(id) {
    // Convert id to an integer
    id = parseInt(id, 10);
    // Retrieve cached atoms
    const atoms = this._cache[id];
    const prefetch_shape = [-5, 15];
    const request_shape = [-20, 20];
    // TODO: adapt reuqest_shape and prefetch_shape based on if a step is requested or it is playing
    // could also first request the frame and then +/-

    // Check if id +/- request_shape is in the cache or requested
    let request = false;
    for (let i = prefetch_shape[0]; i <= prefetch_shape[1]; i++) {
      if (
        this._cache[id + i] === undefined &&
        !this._requested[id + i] &&
        id + i >= 0 &&
        id + i < this._length
      ) {
        request = true;
      }
    }

    // Check if atoms are not cached and the id has not been requested
    if ((atoms === undefined && !this._requested[id]) || request) {
      this._requested[id] = true;
      this.socket.emit("room:frames:get", [id], (data) => {
        this.setFrames(data);
      });
      let new_steps = [];
      for (let i = request_shape[0]; i <= request_shape[1]; i++) {
        new_steps.push(id + i);
      }
      // Filter out negative steps and those already in the cache
      new_steps = new_steps.filter(
        (x) => x >= 0 && this._cache[x] === undefined,
      );
      // Filter out steps larger than the length
      new_steps = new_steps.filter((x) => x < this._length);
      // Log the steps to request
      console.log("Requesting", new_steps);
      // add all requested steps to the requested id list
      new_steps.forEach((step) => (this._requested[step] = true));

      // Emit request for the steps
      // this.socket.emit("room:frames", new_steps, (ack) => {});
      this.socket.emit("room:frames:get", new_steps, (data) => {
        this.setFrames(data);
      });
    }

    return atoms;
  }

  get_length() {
    if (!this._size_timeout) {
      this.socket.emit("room:length:get", (data) => {
        this._length = data;
        const slider = document.getElementById("frameProgress");
        slider.ariaValueMax = Math.max(0, this._length - 1);

        document.getElementById("info").innerHTML =
          `${slider.ariaValueNow} / ${slider.ariaValueMax}`;

        console.log("Cache length is now", this._length);
      });
      this._size_timeout = true;
      setTimeout(() => {
        this._size_timeout = false;
      }, 1000);
    }
    return Math.max(0, this._length);
  }

  getAllAtoms() {
    return this._cache;
  }

  attachWorld(world) {
    this.world = world;
  }

  setFrames(data) {
    const indicesToDelete = [];
    for (const [index, atoms] of Object.entries(data)) {
      if (atoms === null) {
        delete this._cache[index];
      } else {
        this._cache[index] = new Atoms({
          positions: atoms.positions,
          cell: atoms.cell,
          numbers: atoms.numbers,
          colors: atoms.arrays.colors,
          radii: atoms.arrays.radii,
          connectivity: atoms.connectivity,
          calc: atoms.calc,
          pbc: atoms.pbc,
        });
      }
      delete this._requested[index];
    }
    if (!this._received_data) {
      this._received_data = true;
      this.world.setStep(this.world.getStep());
    }
  }
}

export { Cache };
