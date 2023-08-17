// Interface for the communication with Python to retrieve atoms

class Atoms {
    constructor(positions, cell, numbers) {
        this.positions = positions;
        this.cell = cell;
        this.numbers = numbers;
    }
}

class Cache {
    constructor(socket) {
        this._socket = socket;
        this._cache = {};
    }

    async get(id) {
        // if the id is not in the cache, request it from the server
        if (!(id in this._cache)) {
            // create a promise that resolves when the server response is received
            const dataPromise = new Promise((resolve) => {
                this._socket.emit("configuration:id", { id: id }, (data) => {
                    Object.keys(data).forEach((key) => {
                        this._cache[key] = new Atoms(
                            data[key].positions,
                            data[key].cell,
                            data[key].numbers
                        );
                    });
                    resolve();
                });
            });
            // wait for the server response before returning the cached data
            await dataPromise;
        }
        return this._cache[id];
    }

    reset() {
        this._cache = {};
    }
}

export { Cache };
