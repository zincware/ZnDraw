import * as THREE from "three";
import { CSS2DObject } from "three/examples/jsm/renderers/CSS2DRenderer.js";

export const materials = {
  MeshBasicMaterial: new THREE.MeshBasicMaterial({ color: "#ffa500" }),
  MeshLambertMaterial: new THREE.MeshLambertMaterial({ color: "#ffa500" }),
  MeshMatcapMaterial: new THREE.MeshMatcapMaterial({ color: "#ffa500" }),
  MeshPhongMaterial: new THREE.MeshPhongMaterial({
    color: "#ffa500",
    shininess: 100,
  }),
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
  const mesh = new THREE.Mesh(geometry, material);
  const bond = new THREE.Group();
  mesh.userData.type = "bond";
  bond.click = function () {
    this.parent.click();
  };
  bond.add(mesh);
  return bond;
}

function updateBondOrientation(bond, pointX, pointY) {
  const direction = new THREE.Vector3();
  direction.subVectors(pointY, pointX);
  bond.lookAt(pointY);
  const scale =
    direction.length() / 2 / bond.children[0].geometry.parameters.height;
  bond.scale.set(1, 1, scale);
}

function updateParticleScaleAndMaterial(particle, radius, material) {
  const scale = radius / particle.children[0].geometry.parameters.radius;
  particle.children[0].scale.set(scale, scale, scale);
  particle.traverse((x) => {
    if (x instanceof THREE.Mesh) {
      x.material = material;
    }
  });
}

const halfCylinderGeometry = halfCylinderGeometryFactory();
const sphereGeometry = sphereGeometryFactory();
export const speciesMaterial = speciesMaterialFactory();

export function createIndexGroup(particleGroup) {
  const indexGroup = new THREE.Group();

  indexGroup.show = function () {
    if (indexGroup.children.length > 0) {
      return;
    }

    particleGroup.children.forEach((particleSubGroup) => {
      const particle = particleSubGroup.children[0];

      const text = document.createElement("div");
      text.className = "label";
      text.style.color = "black";
      text.textContent = particleSubGroup.name;
      text.style.fontSize = "20px";

      const label = new CSS2DObject(text);
      label.name = `label-${particleSubGroup.name}`;
      label.position.copy(particleSubGroup.position);
      indexGroup.add(label);
    });
  };

  indexGroup.hide = function () {
    while (indexGroup.children.length > 0) {
      indexGroup.remove(indexGroup.children[0]);
    }
  };

  indexGroup.tick = (data) => {
    particleGroup.children.forEach((particleSubGroup) => {
      const label = indexGroup.getObjectByName(
        `label-${particleSubGroup.name}`,
      );
      if (label) {
        label.position.copy(particleSubGroup.position);
      } else {
        indexGroup.remove(label);
      }
    });
  };

  return indexGroup;
}
export function createParticleGroup(config) {
  const particleGroup = new THREE.Group();
  particleGroup.name = "particleGroup";

  particleGroup.get_center = () => {
    const items = [];
    if (config.selected.length > 0) {
      config.selected.forEach((x) => {
        items.push(particleGroup.getObjectByName(x));
      });
    } else {
      particleGroup.children.forEach((x) => {
        items.push(x);
      });
    }

    const center = items
      .reduce((a, b) => a.add(b.position), new THREE.Vector3())
      .divideScalar(items.length);
    return center;
  };

  particleGroup.tick = (data) => {
    if (data == null) {
      return;
    }
    const particles = data["particles"];
    const bonds = data["bonds"];

    // create particle arrays
    let existing_particles = [];
    let new_particles = [];
    let deleted_particles = [];

    particles.forEach((x) => {
      const particleObject = particleGroup.getObjectByName(x.id);

      if (particleObject) {
        existing_particles.push(x);
      } else {
        new_particles.push(x);
      }
    });

    deleted_particles = particleGroup.children.filter(
      (x) => !particles.find((y) => y.id === x.name),
    );

    // create bond arrays
    const all_bonds = particleGroup.children.flatMap((particleSubGroup) =>
      particleSubGroup.children.slice(1),
    );

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
      // TODO this should not depend on FPS!
      if (config.selected.includes(particle.id)) {
        material = speciesMaterial(
          config.config.material,
          "#ffa500",
          config.config.material_wireframe,
        );
      }

      updateParticleScaleAndMaterial(
        particleSubGroup,
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
      particle_mesh.userData.type = "particle";
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
        const [particle1Name, particle2Name] = bond;

        const particle1SubGroup = particleGroup.getObjectByName(particle1Name);
        const particle2SubGroup = particleGroup.getObjectByName(particle2Name);
        const particle1 = particle1SubGroup.children[0];
        const particle2 = particle2SubGroup.children[0];

        const node1 = new THREE.Vector3();
        const node2 = new THREE.Vector3();

        const createBond = (startNode, endNode, startMaterial, name) => {
          const bond_mesh = halfCylinderMesh(
            startNode,
            endNode,
            startMaterial,
            config.config,
          );
          bond_mesh.tick = () => {
            particle1.getWorldPosition(node1);
            particle2.getWorldPosition(node2);
            updateBondOrientation(bond_mesh, startNode, endNode);
          };
          bond_mesh.set_order = (order) => {
            bond_mesh.order = order;
            if (bond_mesh.order == 1) {
              if (bond_mesh.children.length != 1) {
                const bond1 = halfCylinderMesh(
                  startNode,
                  endNode,
                  startMaterial,
                  config.config,
                ).children[0];
                // differentiate between double / triple bond rescale

                // remove all chidren from bond_mesh
                for (let i = bond_mesh.children.length - 1; i >= 0; i--) {
                  bond_mesh.children[i].removeFromParent();
                }

                bond_mesh.add(bond1);
              }
            } else if (bond_mesh.order == 2) {
              if (bond_mesh.children.length != 2) {
                const bond1 = halfCylinderMesh(
                  startNode,
                  endNode,
                  startMaterial,
                  config.config,
                ).children[0];
                const bond2 = bond1.clone();
                bond1.scale.set(0.8, 0.8, 1);
                bond2.scale.set(0.8, 0.8, 1);

                // remove all chidren from bond_mesh
                for (let i = bond_mesh.children.length - 1; i >= 0; i--) {
                  bond_mesh.children[i].removeFromParent();
                }

                bond1.translateX(0.2);
                bond2.translateX(-0.2);
                bond_mesh.add(bond1, bond2);
              }
            } else if (bond_mesh.order == 3) {
              if (bond_mesh.children.length != 3) {
                const bond1 = halfCylinderMesh(
                  startNode,
                  endNode,
                  startMaterial,
                  config.config,
                ).children[0];
                const bond2 = bond1.clone();
                const bond3 = bond2.clone();
                bond1.scale.set(0.7, 0.7, 1);
                bond2.scale.set(0.7, 0.7, 1);
                bond3.scale.set(0.7, 0.7, 1);

                // remove all chidren from bond_mesh
                for (let i = bond_mesh.children.length - 1; i >= 0; i--) {
                  bond_mesh.children[i].removeFromParent();
                }

                bond1.translateX(0.25);
                bond2.translateX(-0.25);
                bond_mesh.add(bond1, bond2, bond3);
              }
            }
          };
          bond_mesh.name = name;
          return bond_mesh;
        };

        const bond_1 = createBond(
          node1,
          node2,
          particle1.material,
          `${particle1Name}-${particle2Name}`,
        );
        const bond_2 = createBond(
          node2,
          node1,
          particle2.material,
          `${particle2Name}-${particle1Name}`,
        );

        particle1SubGroup.add(bond_1);
        particle2SubGroup.add(bond_2);

        bond_1.tick();
        bond_2.tick();
      });

      // update existing bonds
      existing_bonds.forEach((bond) => {
        // get the entry from bonds
        let order = bonds.find(
          (x) =>
            x[0] + "-" + x[1] === bond.name || x[1] + "-" + x[0] === bond.name,
        );
        bond.set_order(order[2]);
        bond.tick();
      });
    }
  };

  return particleGroup;
}
