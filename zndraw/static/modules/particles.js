import * as THREE from 'three';

// each entry in the particleGroup again is a group of [particle, bond1, bond2, bond3, ...]
export const particleGroup = new THREE.Group();


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

const direction = new THREE.Vector3();

// TODO reuse geometry and material for all atoms and just modify the meshes

// a simple memoized function to add something
const halfCylinderGeometryFactory = () => {
    let cache = {};
    let key = ""
    return (bond_size, resolution) => {
        key = bond_size + "_" + resolution;

        if (key in cache) {
            return cache[key];
        }
        else {
            const geometry = new THREE.CylinderGeometry(0.15 * bond_size, 0.15 * bond_size, 1, resolution * 2);
            // shift it so one end rests on the origin
            geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(0, 1 / 2, 0));
            // rotate it the right way for lookAt to work
            geometry.applyMatrix4(new THREE.Matrix4().makeRotationX(THREE.MathUtils.degToRad(90)));
            cache[key] = geometry;
            return geometry;
        }
    }
}

const sphereGeometryFactory = () => {
    let cache = {};
    let key = ""
    return (sphere_size, resolution) => {
        key = sphere_size + "_" + resolution;

        if (key in cache) {
            return cache[key];
        }
        else {
            const geometry = new THREE.SphereGeometry(sphere_size, resolution * 4, resolution * 2);
            cache[key] = geometry;
            return geometry;
        }
    }
}


const speciesMaterialFactory = () => {
    let cache = {};
    let key = ""
    return (name, color, wireframe) => {
        key = name + "_" + color + "_" + wireframe;

        if (key in cache) {
            console.log('Fetching material from cache');
            return cache[key];
        }
        else {
            console.log('Creating new material');
            const material = materials[name].clone()
            material.color.set(color);
            material.wireframe = wireframe;
            cache[key] = material;
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

// TODO rename clean AtomsBonds
// this does not remove the scene!
export function cleanScene(scene) {
    while (particleGroup.children.length > 0) {
        scene.remove(particleGroup.children.shift());
    };
}

export function drawAtoms(atoms, bonds, config, scene) {
    cleanScene(scene);

    atoms.forEach((item) => { addAtom(item, config) });

    if (config["bond_size"] > 0) {
        bonds.forEach((item) => { addBond(item, config) });
    }
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