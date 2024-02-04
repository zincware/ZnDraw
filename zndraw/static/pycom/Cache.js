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

  setFrames(data) {
    const slider = document.getElementById("frameProgress");
    for (const [index, atoms] of Object.entries(data)) {
      if (atoms === undefined) {
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
        slider.ariaValueMax = Object.keys(this._cache).length - 1;
        this.world.setStep(this.world.getStep(), false);
      }
    }
  }
}

export { Cache };
