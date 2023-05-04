import * as THREE from 'three';

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

export function createParticleGroup() {
    const particleGroup = new THREE.Group();

    particleGroup.tick = (data) => {
        if (data == null) {
            return;
        }
        const particles = data["particles"];

        particles.forEach((item) => {
            if (particleGroup.getObjectByName(item.id)) {
                particleGroup.getObjectByName(item.id).position.set(item.x, item.y, item.z);
                // Update size and color if changed
                particleGroup.getObjectByName(item.id).scale.set(item.radius, item.radius, item.radius);
                particleGroup.getObjectByName(item.id).material.color.set(item.color);
            } else {
                const particle = new THREE.Mesh(
                    new THREE.SphereGeometry(1, 32, 32),
                    new THREE.MeshPhongMaterial({ color: item.color })
                );
                particle.position.set(item.x, item.y, item.z);
                particle.scale.set(item.radius, item.radius, item.radius)
                particle.name = item.id;
                particleGroup.add(particle);
            }
        });
        // remove particles that are not in data
        particleGroup.children.forEach((particle) => {
            if (particles.filter((item) => item.id === particle.name).length === 0) {
                particleGroup.remove(particle);
            }
        }
        );
    }

    return particleGroup;
}

function halfCylinderMesh(pointX, pointY, material, config) {
    // // Make a mesh with the geometry
    const direction = new THREE.Vector3();

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

function addBond(particleGroup, item) {
    let particle1 = particleGroup.getObjectByName(item[0]);
    let particle2 = particleGroup.getObjectByName(item[1]);

    let node1 = new THREE.Vector3();
    let node2 = new THREE.Vector3();

    particle1.getWorldPosition(node1);
    particle2.getWorldPosition(node2);

    // temporary config
    let config = {"bond_size": 1.0, "resolution": 8}

    const bond_1 = halfCylinderMesh(node1, node2, particle1.material, config);
    const bond_2 = halfCylinderMesh(node2, node1, particle2.material, config);

    // the atom to look at
    bond_1.name = item[0] + "-" + item[1];
    bond_2.name = item[1] + "-" + item[0];

    // return new THREE.Group().add(bond_1, bond_2);
    return bond_1;
}

export function createBondsGroup(particleGroup) {
    const bondsGroup = new THREE.Group();

    bondsGroup.tick = (data) => {
        if (data == null) {
            return;
        }
        if (bondsGroup.children.length > 10) {
            return; // do not update anymore
        }
        console.log("updating bonds")
        // Update bonds

        const bonds = data["bonds"];

        // remove bonds to particles that no longer exist
        bondsGroup.children.forEach((bond) => {
            if (bonds.filter((item) => item[0] === bond.name.split("-")[0] && item[1] === bond.name.split("-")[1]).length === 0) {
                bondsGroup.remove(bond);
            }
        });

        bonds.forEach((item) => {
            if (bondsGroup.getObjectByName(item[0] + "-" + item[1])) {
                return;
                // console.log("update");
            } else {
                bondsGroup.add(addBond(particleGroup, item));
                // console.log(item);
            }
        });

    }
    return bondsGroup;

}