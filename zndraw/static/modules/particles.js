import * as THREE from 'three';

// each entry in the particleGroup again is a group of [particle, bond1, bond2, bond3, ...]
export const particleGroup = new THREE.Group();
export const arrowGroup = new THREE.Group();


export const materials = {
    "MeshBasicMaterial": new THREE.MeshBasicMaterial({ color: "#ffa500" }),
    "MeshLambertMaterial": new THREE.MeshLambertMaterial({ color: "#ffa500" }),
    "MeshMatcapMaterial": new THREE.MeshMatcapMaterial({ color: "#ffa500" }),
    "MeshPhongMaterial": new THREE.MeshPhongMaterial({ color: "#ffa500" }),
    "MeshPhysicalMaterial": new THREE.MeshPhysicalMaterial({ color: "#ffa500" }),
    "MeshStandardMaterial": new THREE.MeshStandardMaterial({ color: "#ffa500" }),
    "MeshToonMaterial": new THREE.MeshToonMaterial({ color: "#ffa500" }),

};

let node1 = new THREE.Vector3();
let node2 = new THREE.Vector3();
export let box;
export let boxGeometry;

const direction = new THREE.Vector3();

let sphereGeometryFactoryCache = {};
let speciesMaterialFactoryCache = {};
let halfCylinderGeometryFactoryCache = {};

// TODO reuse geometry and material for all atoms and just modify the meshes

// a simple memoized function to add something
const halfCylinderGeometryFactory = () => {
    let key = ""
    return (bond_size, resolution) => {
        key = bond_size + "_" + resolution;

        if (key in halfCylinderGeometryFactoryCache) {
            return halfCylinderGeometryFactoryCache[key];
        }
        else {
            const geometry = new THREE.CylinderGeometry(0.15 * bond_size, 0.15 * bond_size, 1, resolution * 2, 1, true);
            // shift it so one end rests on the origin
            geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(0, 1 / 2, 0));
            // rotate it the right way for lookAt to work
            geometry.applyMatrix4(new THREE.Matrix4().makeRotationX(THREE.MathUtils.degToRad(90)));
            halfCylinderGeometryFactoryCache[key] = geometry;
            return geometry;
        }
    }
}

const sphereGeometryFactory = () => {
    let key = ""
    return (sphere_size, resolution) => {
        key = sphere_size + "_" + resolution;

        if (key in sphereGeometryFactoryCache) {
            return sphereGeometryFactoryCache[key];
        }
        else {
            const geometry = new THREE.SphereGeometry(sphere_size, resolution * 4, resolution * 2);
            sphereGeometryFactoryCache[key] = geometry;
            return geometry;
        }
    }
}


const speciesMaterialFactory = () => {
    let key = ""
    return (name, color, wireframe) => {
        key = name + "_" + color + "_" + wireframe;

        if (key in speciesMaterialFactoryCache) {
            return speciesMaterialFactoryCache[key];
        }
        else {
            console.log('Creating new material');
            const material = materials[name].clone()
            material.color.set(color);
            material.wireframe = wireframe;
            speciesMaterialFactoryCache[key] = material;
            return material;
        }
    }
}

const halfCylinderGeometry = halfCylinderGeometryFactory();
const sphereGeometry = sphereGeometryFactory();
export const speciesMaterial = speciesMaterialFactory();

function halfCylinderMesh(pointX, pointY, material, config) {
    // // Make a mesh with the geometry
    direction.subVectors(pointY, pointX);

    let geometry = halfCylinderGeometry(config["bond_size"], config["resolution"]);

    let mesh = new THREE.Mesh(geometry, material);
    // // Position it where we want
    mesh.position.copy(pointX);
    // // And make it point to where we want
    mesh.lookAt(pointY);
    let scale = (direction.length() / 2) / mesh.geometry.parameters.height;
    mesh.scale.set(1, 1, scale);
    return mesh;
}

function addAtom(item, config) {
    let geometry = sphereGeometry(config["sphere_size"], config["resolution"]);
    let material = speciesMaterial(
        document.getElementById('materialSelect').value, item["color"],
        document.getElementById('wireframe').checked
    );

    const particle = new THREE.Mesh(geometry, material);

    let particle_grp = new THREE.Group();
    particle_grp.add(particle);
    particleGroup.add(particle_grp);

    particle.scale.set(item["radius"], item["radius"], item["radius"]);

    particle.position.set(...item["position"]);
    particle.userData["id"] = item["id"];
    particle.userData["color"] = item["color"];
    particle.userData["bond_ids"] = [];
}

function addBond(item, config) {
    let particle1 = particleGroup.children[item[0]].children[0];
    let particle2 = particleGroup.children[item[1]].children[0];
    particle1.getWorldPosition(node1);
    particle2.getWorldPosition(node2);

    const bond_1 = halfCylinderMesh(node1, node2, particle1.material, config);
    const bond_2 = halfCylinderMesh(node2, node1, particle2.material, config);

    // the atom to look at
    bond_1.userData["target_atom"] = item[1];
    bond_2.userData["target_atom"] = item[0];

    for (let i = 0; i < particleGroup.children.length; i++) {
        if (particleGroup.children[i].children[0].userData["id"] == item[0]) {
            particleGroup.children[i].add(bond_1);
        }
        if (particleGroup.children[i].children[0].userData["id"] == item[1]) {
            particleGroup.children[i].add(bond_2);
        }
    }

}


export function countBonds() {
    let count = 0;
    particleGroup.traverse((child) => {
        if (child.isMesh) {
            count += 1;
        }
    });
    return (count - particleGroup.children.length) / 2; // subtract the number of particles and account for bonds being two meshes
}

function getAtomGrpById(atom_id) {
    let atom = undefined;
    for (let i = 0; i < particleGroup.children.length; i++) {
        if (particleGroup.children[i].children[0].userData["id"] == atom_id) {
            atom = particleGroup.children[i];
            break;
        }
    }
    return atom;
}

/**
 * Get the atom from the ID given by flask
 * @param {BigInt} atom_id 
 * @returns {THREE.Mesh}
 */
export function getAtomById(atom_id) {
    let atom_grp = getAtomGrpById(atom_id);
    if (atom_grp == undefined) {
        return undefined;
    }
    return atom_grp.children[0];
}

/**
 * Compute the center of a list of atoms
 * @param {list} atom_ids 
 */
export function getAtomsCenter(atom_ids) {
    let center = new THREE.Vector3(0, 0, 0);
    let count = 0;
    for (let i = 0; i < atom_ids.length; i++) {
        let atom = getAtomById(atom_ids[i]);
        if (atom == undefined) {
            continue;
        }
        center.add(atom.position);
        count += 1;
    }
    center.divideScalar(count);
    return center;
}

// TODO rename clean AtomsBonds
// this does not remove the scene!
// TODO add setup function that loads the scene inside this module
export function cleanScene(scene) {
    // speciesMaterialFactoryCache = {};
    // sphereGeometryFactoryCache = {};
    // halfCylinderGeometryFactoryCache = {};

    while (particleGroup.children.length > 0) {
        scene.remove(particleGroup.children.shift());
    };

    while (arrowGroup.children.length > 0) {
        scene.remove(arrowGroup.children.shift());
    };
}

export function drawAtoms(atoms, bonds, config, scene) {
    cleanScene(scene);

    atoms.forEach((item) => { addAtom(item, config) });

    if (config["bond_size"] > 0) {
        bonds.forEach((item) => { addBond(item, config) });
    }
    createArrowPerParticle();
    scene.add(particleGroup);
}


/**
 * Update the positions of the particles
 * @param {Float32List} positions as (n_particles, 3)
 */
export function updateParticlePositions(positions) {
    positions.forEach(function (item, index) {
        let per_atom_grp = getAtomGrpById(index);
        let atom = per_atom_grp.children[0];
        node1.set(...item);
        atom.position.copy(node1);

        for (let j = 1; j < per_atom_grp.children.length; j++) {
            let bond = per_atom_grp.children[j];
            node2.set(...positions[bond.userData["target_atom"]]);
            direction.subVectors(node1, node2);
            let scale = (direction.length() / 2);
            if (scale > 1.5) {
                scale = 0.0;
            }
            bond.scale.set(1, 1, scale);
            bond.position.copy(node1);
            bond.lookAt(node2);
        }
    });
};

const raycaster = new THREE.Raycaster();
const pointer = new THREE.Vector2();


/**
 * Print the particle Indices as 2D overlay
 * @param {THREE.Camera} camera 
 * @returns 
 */
export function printIndices(camera) {

    let ids = document.getElementsByClassName("particle-id")
    if (ids.length > 0) {
        return;
    }

    let positions = [];
    let distances = [];
    particleGroup.children.forEach(function (atoms_grp) {
        let item = atoms_grp.children[0];
        let vector = item.position.clone().project(camera);
        vector.x = (vector.x + 1) / 2 * window.innerWidth;
        vector.y = -(vector.y - 1) / 2 * window.innerHeight;
        // if x smaller 0 or larger window width return
        if (vector.x < 50 || vector.x > window.innerWidth - 50) {
            return;
        }
        // if y smaller 0 or larger window height return
        if (vector.y < 50 || vector.y > window.innerHeight - 50) {
            return;
        }


        // between -1 and 1
        pointer.x = (vector.x / window.innerWidth) * 2 - 1;
        pointer.y = - (vector.y / window.innerHeight) * 2 + 1;

        raycaster.setFromCamera(pointer, camera);
        let intersects = raycaster.intersectObjects(particleGroup.children, true);

        if (!(intersects[0].object == item)) {
            return;
        }
        positions.push([vector, item.userData["id"]]);
        distances.push(intersects[0].distance);
    });

    positions.forEach(function (item, index) {

        var text2 = document.createElement('div');
        text2.classList.add("particle-id", "rounded");
        text2.style.position = 'absolute';
        text2.style.fontSize = Math.max(15, parseInt(50 - 0.3 * (distances[index] * Math.max(...distances)))) + 'px';
        text2.style.backgroundColor = "#cccccc";
        text2.innerHTML = item[1];
        text2.style.top = item[0].y + 'px';
        text2.style.left = item[0].x + 'px';
        document.body.appendChild(text2);
    });
}

/**
 * Draw an arrow on each particle with the given direction / magnitude
 * @param {THREE.Vector3} vectors 
 */
export async function updateArrows(vectors){

    arrowGroup.children.forEach((item, index) => {
        const vector = new THREE.Vector3(...vectors[index]);
        item.setDirection(vector.clone().normalize());
        item.setLength(vector.length()); // find a good factor here
        item.position.copy(getAtomById(index).position);
    });
}

export function createArrowPerParticle() {
    particleGroup.children.forEach(function (atoms_grp) {
        let item = atoms_grp.children[0];
        let arrow = new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), item.position, 0, 0xff0000);
        arrowGroup.add(arrow);
    });
    return arrowGroup;
}

export function createBox(size) {
    const geometry = new THREE.BoxGeometry( 1, 1, 1 );
    geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(1/2, 1/2, 1/2));
    boxGeometry = new THREE.EdgesGeometry( geometry );
    const material = new THREE.LineBasicMaterial( { color: 0x000000 } );
    box = new THREE.LineSegments( boxGeometry, material );
    boxGeometry.scale(...size);
    return box;
}

export function updateBox(size) {
    console.log("update box not implemented yet")
}
