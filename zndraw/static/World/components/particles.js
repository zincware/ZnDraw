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

function halfCylinderMesh(pointX, pointY, material, bond_size, resolution) {
  const geometry = halfCylinderGeometry(bond_size, resolution);
  return new THREE.Mesh(geometry, material);
}

const halfCylinderGeometry = halfCylinderGeometryFactory();
const sphereGeometry = sphereGeometryFactory();
export const speciesMaterial = speciesMaterialFactory();

/**
 * Contain a single Particle and its connections
 */
class ParticleGroup extends THREE.Group {
  constructor(particle) {
    super();
    const particle_mesh = new THREE.Mesh(
      sphereGeometry(particle.radius, 10),
      speciesMaterial("MeshPhongMaterial", particle.color, false),
    );
    this.add(particle_mesh);
    this.name = particle.id;
    this._original_material = particle_mesh.material;

    this.position.set(...particle.position);
  }

  update(particle) {
    const scale = particle.radius / this.children[0].geometry.parameters.radius;
    const material = speciesMaterial(
      "MeshPhongMaterial",
      particle.color,
      false,
    );
    // this command may update the selected material
    this._original_material = material;

    this.children[0].scale.set(scale, scale, scale);
    this.position.set(...particle.position);
    this.children.forEach((x) => (x.material = material));
  }

  updateBondOrientation(bond, particle_group) {
    const node1 = new THREE.Vector3();
    const node2 = new THREE.Vector3();

    this.children[0].getWorldPosition(node1);
    particle_group.children[0].getWorldPosition(node2);

    const direction = new THREE.Vector3();
    direction.subVectors(node1, node2);
    bond.lookAt(node2);
    const scale = direction.length() / 2 / bond.geometry.parameters.height;
    bond.scale.set(1, 1, scale);
  }

  connect(particle_group) {
    const bond_mesh = halfCylinderMesh(
      this.children[0],
      particle_group.children[0],
      this.children[0].material,
      1.3,
      8,
    );
    this.add(bond_mesh);
    // Store all the bond information, don't pass the mesh or group here
    this.updateBondOrientation(bond_mesh, particle_group);
  }

  set_selection(selected) {
    if (selected) {
      // see https://threejs.org/examples/#webgl_postprocessing_unreal_bloom_selective for inspiration
      const material = speciesMaterial("MeshPhongMaterial", "#ffa500", false);
      this.children.forEach((x) => (x.material = material));
    } else {
      this.children.forEach((x) => (x.material = this._original_material));
    }
  }
}

/**
 * Contain all Particles of the World.
 */
class ParticlesGroup extends THREE.Group {
  constructor(socket, cache) {
    super();
    this.name = "particlesGroup";
    this.cache = cache;
    socket.on("config", (data) => {
      this.config = data;
      console.log("ParticlesGroup config:");
      console.log(this.config);
    });
  }

  tick() {}

  _updateParticles(particles) {
    let existing_particles = [];
    let new_particles = [];
    let deleted_particles = [];

    for (const x of particles) {
      const particleObject = this.getObjectByName(x.id);

      if (particleObject) {
        existing_particles.push(x);
      } else {
        new_particles.push(x);
      }
    }

    // get all particles who have an id larger than particles.length
    deleted_particles = this.children.filter((x) => x.name >= particles.length);

    new_particles.forEach((particle) => {
      this.add(new ParticleGroup(particle));
    });

    existing_particles.forEach((particle) => {
      const particleSubGroup = this.getObjectByName(particle.id);
      particleSubGroup.update(particle);
    });

    // remove deleted particles
    deleted_particles.forEach((particle) => {
      particle.removeFromParent();
      // config.selected = config.selected.filter((e) => e !== particle.name);
    });

    // console.log("Having existing particles: " + existing_particles.length + " and adding " + new_particles.length + " and removing " + deleted_particles.length);
  }

  _updateBonds(bonds) {
    // create bond arrays
    const all_bonds = this.children.flatMap((particleSubGroup) =>
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

    deleted_bonds.forEach((bond) => {
      bond.removeFromParent();
    });

    new_bonds.forEach((bond) => {
      const [particle1Name, particle2Name] = bond;
      const particle1SubGroup = this.getObjectByName(particle1Name);
      const particle2SubGroup = this.getObjectByName(particle2Name);

      particle1SubGroup.connect(particle2SubGroup);
      particle2SubGroup.connect(particle1SubGroup);
    });
  }

  step(frame) {
    const particles = this.cache.get(frame);
    if (particles == null) {
      // nothing to display
    } else {
      this._updateParticles(particles);
      this._updateBonds(particles.connectivity);
    }
  }
}

export { ParticlesGroup };
