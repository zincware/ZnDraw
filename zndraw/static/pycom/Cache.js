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
    constructor({ positions, cell, numbers, colors, radii }) {
        this.positions = positions;
        this.cell = cell;
        this.numbers = numbers;
        this.colors = colors;
        this.radii = radii;
    }

    [Symbol.iterator]() {
        let index = 0;
        const atoms = this; // this does not work
        return {
            next() {
                if (index < atoms.positions.length) {
                    const atom = new Atom(
                        {
                            position: atoms.positions[index],
                            number: atoms.numbers[index],
                            color: atoms.colors[index],
                            radius: atoms.radii[index],
                            id: index
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
                            radii: data[key].radii
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
        const value = this._cache[id];
        if (value === undefined) {
            this._get(id);
        }
        return value;
    }

    reset() {
        this._cache = {};
    }
}

export { Cache };
