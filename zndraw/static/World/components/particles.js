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
      // console.log("Creating new material");
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

    // create particle arrays
    let existing_particles = particles.filter((x) =>
      particleGroup.getObjectByName(x.id),
    );
    let new_particles = particles.filter(
      (x) => !particleGroup.getObjectByName(x.id),
    );
    let deleted_particles = particleGroup.children.filter(
      (x) => !particles.find((y) => y.id === x.name),
    );

    // create bond arrays
    let all_bonds = [];
    particleGroup.children.forEach((particleSubGroup) => {
      particleSubGroup.children.forEach((item, idx) => {
        if (idx === 0) {
          return;
        }
        all_bonds.push(item);
      });
    });

    let existing_bonds = all_bonds.filter(
      (x) =>
        bonds.find((y) => y[0] + "-" + y[1] === x.name) ||
        bonds.find((y) => y[1] + "-" + y[0] === x.name),
    );
    let new_bonds = bonds.filter(
      (x) =>
        !all_bonds.find((y) => y.name === x[0] + "-" + x[1]) &&
        !all_bonds.find((y) => y.name === x[1] + "-" + x[0]),
    );
    let deleted_bonds = all_bonds.filter(
      (x) =>
        !bonds.find((y) => y[0] + "-" + y[1] === x.name) &&
        !bonds.find((y) => y[1] + "-" + y[0] === x.name),
    );

    // console.log("Having existing particles: " + existing_particles.length + " and adding " + new_particles.length + " and removing " + deleted_particles.length);
    // console.log("Having existing bonds: " + existing_bonds.length + " and adding " + new_bonds.length + " and removing " + deleted_bonds.length);

    // update existing particles
    existing_particles.forEach((particle) => {
      const particleSubGroup = particleGroup.getObjectByName(particle.id);

      particleSubGroup.position.set(particle.x, particle.y, particle.z);

      let material = speciesMaterial(
        config.config.material,
        particle.color,
        config.config.material_wireframe,
      );

      // handle selected particles
      if (config.selected.includes(particle.id)) {
        material = speciesMaterial(
          config.config.material,
          "ffa500",
          config.config.material_wireframe,
        );
      }

      updateParticleScaleAndMaterial(
        particleSubGroup.children[0],
        particle.radius * config.config.sphere_size,
        material,
      );
    });

    // create new particles
    new_particles.forEach((particle) => {
      const particle_mesh = new THREE.Mesh(
        sphereGeometry(particle.radius, config.config.resolution),
        speciesMaterial(
          config.config.material,
          particle.color,
          config.config.material_wireframe,
        ),
      );
      const particleSubGroup = new THREE.Group();
      particleSubGroup.add(particle_mesh);
      particleSubGroup.name = particle.id;
      particleSubGroup.position.set(particle.x, particle.y, particle.z);

      // CLICK EVENT
      particleSubGroup.click = () => {
        if (config.selected.includes(particle.id)) {
          if (!config.pressed_keys.Shift) {
            config.selected = [];
          } else {
            config.selected = config.selected.filter((e) => e !== particle.id);
          }
        } else {
          if (!config.pressed_keys.Shift) {
            config.selected = [particle.id];
          } else {
            config.selected.push(particle.id);
          }
        }
      };

      particleGroup.add(particleSubGroup);
    });

    // remove deleted particles
    deleted_particles.forEach((particle) => {
      particle.removeFromParent();
      config.selected = config.selected.filter((e) => e !== particle.name);
    });

    if (config.config.bond_size > 0) {
      // remove bonds that are not in data
      deleted_bonds.forEach((bond) => {
        bond.removeFromParent();
      });
      new_bonds.forEach((bond) => {
        const particle1SubGroup = particleGroup.getObjectByName(bond[0]);
        const particle2SubGroup = particleGroup.getObjectByName(bond[1]);
        const particle1 = particle1SubGroup.children[0];
        const particle2 = particle2SubGroup.children[0];

        const node1 = new THREE.Vector3();
        const node2 = new THREE.Vector3();

        // const direction = new THREE.Vector3();
        // direction.subVectors(node1, node2);

        const bond_1 = halfCylinderMesh(
          node1,
          node2,
          particle1.material,
          config.config,
        );
        const bond_2 = halfCylinderMesh(
          node2,
          node1,
          particle2.material,
          config.config,
        );

        bond_1.tick = () => {
          particle1.getWorldPosition(node1);
          particle2.getWorldPosition(node2);
          updateBondOrientation(bond_1, node1, node2);
          bond_1.material = particle1.material;
        };
        bond_2.tick = () => {
          particle1.getWorldPosition(node1);
          particle2.getWorldPosition(node2);
          updateBondOrientation(bond_2, node2, node1);
          bond_2.material = particle2.material;
        };

        // the atom to look at
        bond_1.name = bond[0] + "-" + bond[1];
        bond_2.name = bond[1] + "-" + bond[0];

        particle1SubGroup.add(bond_1);
        particle2SubGroup.add(bond_2);

        bond_1.tick();
        bond_2.tick();
      });

      // update existing bonds
      existing_bonds.forEach((bond) => {
        bond.tick();
      });
    }
  };

  return particleGroup;
}
