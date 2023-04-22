import * as THREE from 'three';

export const atomsGroup = new THREE.Group();

export const bondsGroup = new THREE.Group();
export const bondsGroup_1 = new THREE.Group();
export const bondsGroup_2 = new THREE.Group();

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
            console.log('Fetching CylinderGeometry from cache');
            return cache[key];
        }
        else {
            console.log('Creating new CylinderGeometry');
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
            console.log('Fetching SphereGeometry from cache');
            return cache[key];
        }
        else {
            console.log('Creating new SphereGeometry');
            const geometry = new THREE.SphereGeometry(sphere_size, resolution * 4, resolution * 2);
            cache[key] = geometry;
            return geometry;
        }
    }
}

const halfCylinderGeometry = halfCylinderGeometryFactory();
const sphereGeometry = sphereGeometryFactory();

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

    const particle = new THREE.Mesh(geometry, materials[document.getElementById('materialSelect').value].clone());
    atomsGroup.add(particle);
    particle.material.color.set(item["color"]);
    particle.material.wireframe = document.getElementById('wireframe').checked;
    particle.scale.set(item["radius"], item["radius"], item["radius"]);

    particle.position.set(...item["position"]);
    particle.userData["id"] = item["id"];
    particle.userData["color"] = item["color"];
    particle.userData["bond_ids"] = [];
}

function addBond(item, config) {
    atomsGroup.children[item[0]].getWorldPosition(node1);
    atomsGroup.children[item[1]].getWorldPosition(node2);

    const bond_1 = halfCylinderMesh(node1, node2, atomsGroup.children[item[0]].material, config);
    const bond_2 = halfCylinderMesh(node2, node1, atomsGroup.children[item[1]].material, config);

    bond_1.userData["atom_id"] = item[0];
    bond_2.userData["atom_id"] = item[1];

    atomsGroup.children[item[0]].userData["bond_ids"].push(bond_1.id);
    atomsGroup.children[item[0]].userData["bond_ids"].push(bond_2.id);
    atomsGroup.children[item[1]].userData["bond_ids"].push(bond_1.id);
    atomsGroup.children[item[1]].userData["bond_ids"].push(bond_2.id);


    bondsGroup_1.add(bond_1);
    bondsGroup_2.add(bond_2);

}

// TODO rename clean AtomsBonds
// this does not remove the scene!
export function cleanScene(scene) {
    while (atomsGroup.children.length > 0) {
        scene.remove(atomsGroup.children.shift());
    };
    while (bondsGroup_1.children.length > 0) {
        scene.remove(bondsGroup_1.children.shift());
    };
    while (bondsGroup_2.children.length > 0) {
        scene.remove(bondsGroup_2.children.shift());
    };
}

export function drawAtoms(atoms, bonds, config, scene) {
    cleanScene(scene);

    atoms.forEach(function (item, index) {
        // console.log("Adding item " + index + " to scene(" + item + ")");
        addAtom(item, config);
    });

    scene.add(atomsGroup);
    if (config["bond_size"] > 0) {

        bonds.forEach(function (item, index) {
            // console.log("Adding item " + index + " to scene(" + item + ")");
            addBond(item, config);

        });

        bondsGroup.add(bondsGroup_1);
        bondsGroup.add(bondsGroup_2);
        scene.add(bondsGroup);

    }

}