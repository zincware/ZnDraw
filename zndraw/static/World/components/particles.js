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
    key = `${bond_size}_${resolution}`;

    if (key in halfCylinderGeometryFactoryCache) {
      return halfCylinderGeometryFactoryCache[key];
    }
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
  };
};

const sphereGeometryFactory = () => {
  let key = "";
  return (sphere_size, resolution) => {
    key = `${sphere_size}_${resolution}`;

    if (key in sphereGeometryFactoryCache) {
      return sphereGeometryFactoryCache[key];
    }
    const geometry = new THREE.SphereGeometry(
      sphere_size,
      resolution * 4,
      resolution * 2,
    );
    sphereGeometryFactoryCache[key] = geometry;
    return geometry;
  };
};

const speciesMaterialFactory = () => (name, color, wireframe) => {
  const key = `${name}_${color}_${wireframe}`;
  // const key = (name, color, wireframe);

  if (key in speciesMaterialFactoryCache) {
    return speciesMaterialFactoryCache[key];
  }
  // console.log("Creating new material");
  const material = materials[name].clone();
  material.color.set(color);
  material.wireframe = wireframe;
  speciesMaterialFactoryCache[key] = material;
  return material;
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
  constructor(particle, resolution, material, wireframe) {
    super();

    this.resolution = resolution;
    this.material = material;
    this.wireframe = wireframe;

    this.bonds = [];

    const particle_mesh = new THREE.Mesh(
      sphereGeometry(particle.radius, this.resolution),
      speciesMaterial(this.material, particle.color, this.wireframe),
    );
    this.add(particle_mesh);
    this.name = particle.id;
    this._original_material = particle_mesh.material;

    this.position.set(...particle.position);
  }

  update(particle) {
    const scale = particle.radius / this.children[0].geometry.parameters.radius;
    const material = speciesMaterial(
      this.material,
      particle.color,
      this.wireframe,
    );
    // this command may update the selected material
    this._original_material = material;

    this.children[0].scale.set(scale, scale, scale);
    this.position.set(...particle.position);
    this.children.forEach((x) => (x.material = material));
  }

  updateBonds(bonds) {
    const bondsToRemove = [];

    bonds.forEach(([_, targetParticleGroup, bond_type]) => {
      const bond_mesh = this.bonds.find(
        (bond) => bond.particle_group === targetParticleGroup,
      );
      if (bond_mesh) {
        this.updateBondOrientation(bond_mesh.bond, targetParticleGroup);
      } else {
        this.connect(targetParticleGroup);
      }
    });

    this.bonds = this.bonds.filter((bond) => {
      const exists = bonds.some(
        ([_, particle_group]) => particle_group === bond.particle_group,
      );
      if (!exists) {
        bondsToRemove.push(bond);
        return false;
      }
      return true;
    });

    bondsToRemove.forEach((bondToRemove) => {
      this.remove(bondToRemove.bond);
    });
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
      this.resolution,
    );
    this.bonds.push({ bond: bond_mesh, particle_group });
    this.add(bond_mesh);
    // Store all the bond information, don't pass the mesh or group here
    this.updateBondOrientation(bond_mesh, particle_group);
  }

  set_selection(selected) {
    if (selected) {
      // see https://threejs.org/examples/#webgl_postprocessing_unreal_bloom_selective for inspiration
      const material = speciesMaterial(
        this.material,
        "#ffa500",
        this.wireframe,
      );
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
    this.resolution = 10;
    this.material = "MeshPhongMaterial";
    this.wireframe = false;
    this.show_bonds = true;

    this.bonds_exist = false;
  }

  rebuild(resolution, material, wireframe, bonds) {
    // remove all children
    // this.children.forEach((x) => x.removeFromParent());
    this.clear();
    this.resolution = resolution;
    this.material = material;
    this.wireframe = wireframe;
    this.show_bonds = bonds;
  }

  tick() {}

  _updateParticles(particles) {
    const existing_particles = [];
    const new_particles = [];
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
      this.add(
        new ParticleGroup(
          particle,
          this.resolution,
          this.material,
          this.wireframe,
        ),
      );
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
    const getBondsForParticle = (particleName) => {
      const bondsWithGroups = bonds
        .filter(([a, b]) => a === particleName || b === particleName)
        .map(([a, b, bond_type]) =>
          a === particleName
            ? [a, this.getObjectByName(b), bond_type]
            : [b, this.getObjectByName(a), bond_type],
        );
      return bondsWithGroups;
    };

    this.children.forEach((particleSubGroup) => {
      const particleName = particleSubGroup.name;
      const bondsForParticle = getBondsForParticle(particleName);
      particleSubGroup.updateBonds(bondsForParticle);
    });
  }

  step(frame) {
    const particles = this.cache.get(frame);
    if (particles == null) {
      // nothing to display
    } else {
      this._updateParticles(particles);
      if (this.show_bonds) {
        this._updateBonds(particles.connectivity);
      }
    }
  }
}

class CellGroup extends THREE.Group {
  constructor(cache) {
    super();
    this.name = "cellGroup";
    this.cache = cache;
    this.is_visible = false;
  }

  step(frame) {
    if (!this.is_visible) {
      this.clear();
      return;
    }
    const particles = this.cache.get(frame);
    const cell = particles.cell;

    this.clear();
    const boxGeometry = new THREE.BoxGeometry(
      cell[0][0],
      cell[1][1],
      cell[2][2],
    );
    const wireframe = new THREE.EdgesGeometry(boxGeometry);
    this.cell_lines = new THREE.LineSegments(
      wireframe,
      new THREE.LineBasicMaterial({ color: "#000000" }),
    );
    this.cell_lines.position.set(
      cell[0][0] / 2,
      cell[1][1] / 2,
      cell[2][2] / 2,
    );
    // this.cell_lines.set_selection = (selected) => { };
    this.add(this.cell_lines);
  }

  set_visibility(visible) {
    this.is_visible = visible;
  }
}

class ParticleIndexGroup extends THREE.Group {
  constructor(particlesGroup, camera) {
    super();
    this.particlesGroup = particlesGroup;
    this.camera = camera;

    this.show_labels = false;
    this.label_offset = 0;

    document.addEventListener("keydown", (event) => {
      if (event.repeat) {
        return;
      }
      if (event.isComposing || event.key === "i") {
        this.toggle();
      }
    });
  }

  rebuild(label_offset) {
    this.label_offset = label_offset;
    this.clear();
  }

  toggle() {
    this.show_labels = !this.show_labels;
    this.clear();
  }

  step(frame) {
    if (this.show_labels) {
      this.clear();
    }
  }

  tick() {
    if (this.show_labels) {
      const raycaster = new THREE.Raycaster();
      raycaster.camera = this.camera;

      const direction = new THREE.Vector3();

      function get2dPositions(position, camera) {
        const vector = position.clone();
        vector.project(camera);

        const vector2d = new THREE.Vector2(vector.x, vector.y);

        // add shifts +/- x/y if you want to sample from the top/bottom/left/right of the particle

        return [vector2d];
      }

      this.particlesGroup.children.forEach((object) => {
        // combine all intersects from the center and top/bottom/left/right of the particle
        let visible = true;
        let intersects;

        // center
        const positions = get2dPositions(object.position, this.camera);

        positions.forEach((position) => {
          raycaster.setFromCamera(position, this.camera);
          intersects = raycaster.intersectObjects(this.particlesGroup.children);

          if (intersects.length > 0 && intersects[0].object.parent !== object) {
            visible = false;
          }
        });

        if (!visible) {
          // get the `${particle.name}-label`; object from this and remove it
          const label = this.getObjectByName(`${object.name}-label`);
          this.remove(label);
        } else if (!this.getObjectByName(`${object.name}-label`)) {
          // create a div with unicode f00d
          const text = document.createElement("div");
          const blank = "\u00A0";
          const line = "\u23AF";
          text.className = "label";
          if (this.label_offset === 0) {
            text.textContent = object.name;
          } else if (this.label_offset > 0) {
            text.textContent = `${blank.repeat(
              this.label_offset * 2 + 6,
            )}\u00D7${line.repeat(this.label_offset)}${blank}${object.name}`;
          } else {
            text.textContent = `${object.name}${blank}${line.repeat(
              this.label_offset * -1,
            )}\u00D7${blank.repeat(this.label_offset * -2 + 6)}`;
          }
          // textContent = object.name + label_offset * \u23AF + label_offset * \u00A0
          // text.textContent = `\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00D7\u23AF\u23AF\u23AF\u00A0${object.name}`;
          // text.textContent = object.name;
          // text.textContent = String.fromCodePoint(0xf00d);
          text.style.fontSize = "20px";
          // text.style.color = `#${object.children[0].material.color.getHexString()}`;
          text.style.textShadow = "1px 1px 1px #000000";

          // const text = document.createElement('div');
          // // <span class="badge bg-secondary">New</span>
          // text.className = 'label';
          // text.textContent =  `\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0\u00A0X------${object.name}`;
          // text.style.fontSize = '20px';
          // text.style.color = `#${object.children[0].material.color.getHexString()}`
          // // text-shadow: #FC0 1px 0 10px;
          //
          // console.log(object.children[0].material.color);

          const label = new CSS2DObject(text);
          label.position.set(...object.position);
          label.name = `${object.name}-label`;
          this.add(label);
        }
      });
    }
  }
}

export { ParticlesGroup, ParticleIndexGroup, CellGroup };
