import * as THREE from "three";

export const materials = {
  MeshBasicMaterial: new THREE.MeshBasicMaterial({ color: "#ffa500" }),
  MeshLambertMaterial: new THREE.MeshLambertMaterial({ color: "#ffa500" }),
  MeshMatcapMaterial: new THREE.MeshMatcapMaterial({ color: "#ffa500" }),
  MeshPhongMaterial: new THREE.MeshPhongMaterial({ color: "#ffa500" }),
  MeshPhysicalMaterial: new THREE.MeshPhysicalMaterial({ color: "#ffa500" }),
  MeshStandardMaterial: new THREE.MeshStandardMaterial({ color: "#ffa500" }),
  MeshToonMaterial: new THREE.MeshToonMaterial({ color: "#ffa500" }),
};

const sphereGeometryFactoryCache = {};
const speciesMaterialFactoryCache = {};
const halfCylinderGeometryFactoryCache = {};

// TODO reuse geometry and material for all atoms and just modify the meshes

// a simple memoized function to add something
const halfCylinderGeometryFactory = () => {
  let key = "";
  return (bond_size, resolution) => {
    key = bond_size + "_" + resolution;

    if (key in halfCylinderGeometryFactoryCache) {
      return halfCylinderGeometryFactoryCache[key];
    } else {
      const geometry = new THREE.CylinderGeometry(
        0.15 * bond_size,
        0.15 * bond_size,
        1,
        resolution * 2,
        1,
        true,
      );
      // shift it so one end rests on the origin
      geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(0, 1 / 2, 0));
      // rotate it the right way for lookAt to work
      geometry.applyMatrix4(
        new THREE.Matrix4().makeRotationX(THREE.MathUtils.degToRad(90)),
      );
      halfCylinderGeometryFactoryCache[key] = geometry;
      return geometry;
    }
  };
};

const sphereGeometryFactory = () => {
  let key = "";
  return (sphere_size, resolution) => {
    key = sphere_size + "_" + resolution;

    if (key in sphereGeometryFactoryCache) {
      return sphereGeometryFactoryCache[key];
    } else {
      const geometry = new THREE.SphereGeometry(
        sphere_size,
        resolution * 4,
        resolution * 2,
      );
      sphereGeometryFactoryCache[key] = geometry;
      return geometry;
    }
  };
};

const speciesMaterialFactory = () => {
  return (name, color, wireframe) => {
    const key = name + "_" + color + "_" + wireframe;
    // const key = (name, color, wireframe);

    if (key in speciesMaterialFactoryCache) {
      return speciesMaterialFactoryCache[key];
    } else {
      console.log("Creating new material");
      const material = materials[name].clone();
      material.color.set(color);
      material.wireframe = wireframe;
      speciesMaterialFactoryCache[key] = material;
      return material;
    }
  };
};

function halfCylinderMesh(pointX, pointY, material, config) {
  const geometry = halfCylinderGeometry(
    config["bond_size"],
    config["resolution"],
  );
  return new THREE.Mesh(geometry, material);
}

function updateBondOrientation(bond, pointX, pointY) {
  const direction = new THREE.Vector3();
  direction.subVectors(pointY, pointX);
  bond.lookAt(pointY);
  const scale = direction.length() / 2 / bond.geometry.parameters.height;
  bond.scale.set(1, 1, scale);
}

function updateParticleScaleAndMaterial(particle, radius, material) {
  const scale = radius / particle.geometry.parameters.radius;
  particle.scale.set(scale, scale, scale);
  particle.material = material;
}

const halfCylinderGeometry = halfCylinderGeometryFactory();
const sphereGeometry = sphereGeometryFactory();
export const speciesMaterial = speciesMaterialFactory();

export function createParticleGroup(config) {
  const particleGroup = new THREE.Group();

  particleGroup.tick = (data) => {
    if (data == null) {
      return;
    }
    const particles = data["particles"];
    const bonds = data["bonds"];

    // remove bonds that are not in data
    particleGroup.children.forEach((particle) => {
      particle.children.forEach((bond, idx) => {
        if (idx === 0) {
          // entry 0 is the particle itself
          return;
        }
        const bond_a = bonds.filter(
          (item) => item[0] + "-" + item[1] === bond.name,
        );
        const bond_b = bonds.filter(
          (item) => item[1] + "-" + item[0] === bond.name,
        );

        if (bond_a.length === 0 && bond_b.length === 0) {
          console.log("Removing bond " + bond.name + " in " + bonds[0]);
          particle.remove(bond);
        }
      });
    });

    // remove particles that are not in data
    particleGroup.children.forEach((particle) => {
      if (particles.filter((item) => item.id === particle.name).length === 0) {
        particleGroup.remove(particle);
        // remove from config.selected
        config.selected = config.selected.filter(
          (item) => item !== particle.name,
        );
      }
    });
    particles.forEach((item) => {
      if (particleGroup.getObjectByName(item.id)) {
        particleGroup
          .getObjectByName(item.id)
          .position.set(item.x, item.y, item.z);
        // Update size and color if changed
        let material = speciesMaterial(
          config.config.material,
          item.color,
          config.config.material_wireframe,
        );

        // handle selected particles
        if (config.selected.includes(item.id)) {
          material = speciesMaterial(
            config.config.material,
            "ffa500",
            config.config.material_wireframe,
          );
        }

        updateParticleScaleAndMaterial(
          particleGroup.getObjectByName(item.id).children[0],
          item.radius * config.config.sphere_size,
          material,
        );
      } else {
        const particle = new THREE.Mesh(
          sphereGeometry(item.radius, config.config.resolution),
          speciesMaterial(
            config.config.material,
            item.color,
            config.config.material_wireframe,
          ),
        );
        const particleSubGroup = new THREE.Group();
        particleSubGroup.add(particle);
        particleSubGroup.name = item.id;
        particleSubGroup.position.set(item.x, item.y, item.z);

        // CLICK EVENT
        particleSubGroup.click = () => {
          if (config.selected.includes(item.id)) {
            if (!config.pressed_keys.Shift) {
              config.selected = [];
            } else {
              config.selected = config.selected.filter((e) => e !== item.id);
            }
          } else {
            if (!config.pressed_keys.Shift) {
              config.selected = [item.id];
            } else {
              config.selected.push(item.id);
            }
          }
        };

        particleGroup.add(particleSubGroup);
      }
    });

    // now update bonds
    if (config.config.bond_size > 0) {
      bonds.forEach((item) => {
        const particle1SubGroup = particleGroup.getObjectByName(item[0]);
        const particle2SubGroup = particleGroup.getObjectByName(item[1]);
        const particle1 = particle1SubGroup.children[0];
        const particle2 = particle2SubGroup.children[0];

        const node1 = new THREE.Vector3();
        const node2 = new THREE.Vector3();

        particle1.getWorldPosition(node1);
        particle2.getWorldPosition(node2);

        // existing bonds
        let bond_1 = particleGroup.getObjectByName(item[0] + "-" + item[1]);
        let bond_2 = particleGroup.getObjectByName(item[1] + "-" + item[0]);

        const direction = new THREE.Vector3();
        direction.subVectors(node1, node2);

        if (bond_1 && bond_2) {
          // update bond orientation
          updateBondOrientation(bond_1, node1, node2);
          updateBondOrientation(bond_2, node2, node1);
          bond_1.material = particle1.material;
          bond_2.material = particle2.material;
        } else {
          // temporary config
          console.log("Creating new bond");

          bond_1 = halfCylinderMesh(
            node1,
            node2,
            particle1.material,
            config.config,
          );
          bond_2 = halfCylinderMesh(
            node2,
            node1,
            particle2.material,
            config.config,
          );

          // the atom to look at
          bond_1.name = item[0] + "-" + item[1];
          bond_2.name = item[1] + "-" + item[0];

          particle1SubGroup.add(bond_1);
          particle2SubGroup.add(bond_2);

          // update bond orientation
          updateBondOrientation(bond_1, node1, node2);
          updateBondOrientation(bond_2, node2, node1);
        }
      });
    }
  };

  return particleGroup;
}
