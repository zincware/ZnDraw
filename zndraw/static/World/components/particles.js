import * as THREE from "three";
import { CSS2DObject } from "three/examples/jsm/renderers/CSS2DRenderer.js";

const materials = {
  MeshBasicMaterial: new THREE.MeshBasicMaterial(),
  MeshLambertMaterial: new THREE.MeshLambertMaterial(),
  MeshMatcapMaterial: new THREE.MeshMatcapMaterial(),
  MeshPhongMaterial: new THREE.MeshPhongMaterial({
    shininess: 100,
  }),
  MeshPhysicalMaterial: new THREE.MeshPhysicalMaterial(),
  MeshStandardMaterial: new THREE.MeshStandardMaterial(),
  MeshToonMaterial: new THREE.MeshToonMaterial(),
};

/**
 * Contain all Particles of the World.
 */
class ParticlesGroup extends THREE.Group {
  constructor(socket, cache) {
    super();
    this.name = "particlesGroup";
    this.cache = cache;
    this.socket = socket;

    this.bonds_exist = false;
    this.selection = [];

    this.particle_cache = undefined;
    this.bonds_mesh = undefined;

    // rebuild attributes
    this.resolution = 10;
    this.material = "MeshPhongMaterial";
    this.wireframe = false;
    this.show_bonds = true;
    this.particle_size = 1;
    this.bonds_size = 1;
  }

  rebuild(resolution, material, wireframe, bonds, particle_size, bonds_size) {
    this.resolution = resolution;
    this.material = material;
    this.wireframe = wireframe;
    this.show_bonds = bonds;
    this.particle_size = particle_size;
    this.bonds_size = bonds_size;
    this.clear();
    this.particles_mesh = undefined;
    this.bonds_mesh = undefined;
  }

  tick() {}

  click(instanceId, shift, object) {
    if (instanceId !== undefined) {
      if (object === this.bonds_mesh) {
        const all_bonds = this._get_all_bonds(this.particle_cache.connectivity);
        instanceId = all_bonds[instanceId][0];
      }
      if (shift) {
        if (this.selection.includes(instanceId)) {
          this.selection = this.selection.filter((x) => x !== instanceId);
        } else {
          this.selection.push(instanceId);
        }
      } else {
        if (this.selection.includes(instanceId)) {
          this.selection = [];
        } else {
          this.selection = [instanceId];
        }
      }
    } else {
      this.selection = [];
    }
    this.step();
    // trigger this._updateParticles() to update the selection
  }

  _get_particle_mesh(particle) {
    const particles_geometry = new THREE.SphereGeometry(
      this.particle_size,
      this.resolution * 4,
      this.resolution * 2,
    );
    const particles_material = materials[this.material].clone();
    particles_material.wireframe = this.wireframe;

    const particles_mesh = new THREE.InstancedMesh(
      particles_geometry,
      particles_material,
      particle.positions.length,
    );
    return particles_mesh;
  }

  _updateParticles(particles) {
    if (particles === undefined) {
      return;
    }
    if (particles.positions.length === 0) {
      this.remove(this.particles_mesh);
      this.particles_mesh = undefined;
      return;
    }
    if (this.particles_mesh === undefined) {
      this.particles_mesh = this._get_particle_mesh(particles);
      this.add(this.particles_mesh);
    }
    if (this.particles_mesh.count !== particles.positions.length) {
      this.remove(this.particles_mesh);
      this.particles_mesh = this._get_particle_mesh(particles);
      this.add(this.particles_mesh);
    }
    const matrix = new THREE.Matrix4();
    const dummy = new THREE.Object3D();
    const color = new THREE.Color();

    for (let i = 0; i < this.particles_mesh.count; i++) {
      this.particles_mesh.getMatrixAt(i, matrix);
      matrix.decompose(dummy.position, dummy.rotation, dummy.scale);

      dummy.position.set(...particles.positions[i]);
      dummy.scale.x = dummy.scale.y = dummy.scale.z = particles.radii[i];

      dummy.updateMatrix();
      this.particles_mesh.setMatrixAt(i, dummy.matrix);
      if (this.selection.includes(i)) {
        color.set(0xffa500);
      } else {
        color.set(particles.colors[i]);
      }
      this.particles_mesh.setColorAt(i, color);
    }
    this.particles_mesh.instanceMatrix.needsUpdate = true;
    this.particles_mesh.instanceColor.needsUpdate = true;
  }

  _get_bonds_mesh(bonds) {
    const geometry = new THREE.CylinderGeometry(
      0.15 * this.bonds_size,
      0.15 * this.bonds_size,
      1,
      this.resolution * 2,
      1,
      true,
    );
    // shift it so one end rests on the origin
    geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(0, 1 / 2, 0));
    // rotate it the right way for lookAt to work
    geometry.applyMatrix4(
      new THREE.Matrix4().makeRotationX(THREE.MathUtils.degToRad(90)),
    );

    const material = materials[this.material].clone();
    material.wireframe = this.wireframe;
    // bonds A-B from 0..n, bonds B-A from n+1..2n
    const mesh = new THREE.InstancedMesh(geometry, material, bonds.length * 2);
    return mesh;
  }

  _get_all_bonds(bonds) {
    const getBondsForParticle = (instanceId) => {
      return bonds
        .filter(([a, b]) => a === instanceId || b === instanceId)
        .map(([a, b, bond_type]) =>
          a === instanceId ? [a, b, bond_type] : [b, a, bond_type],
        );
    };
    const allBonds = [];
    for (
      let instanceId = 0;
      instanceId < this.particles_mesh.count;
      instanceId++
    ) {
      const bondsForParticle = getBondsForParticle(instanceId);
      bondsForParticle.forEach(([_, targetInstanceId, bond_type]) => {
        allBonds.push([instanceId, targetInstanceId, bond_type]);
      });
    }
    return allBonds;
  }

  _updateBonds(bonds) {
    if (bonds === undefined) {
      return;
    }
    if (bonds.length === 0) {
      this.remove(this.bonds_mesh);
      this.bonds_mesh = undefined;
      return;
    }
    if (this.bonds_mesh === undefined) {
      this.bonds_mesh = this._get_bonds_mesh(bonds);
      this.add(this.bonds_mesh);
    }
    if (this.bonds_mesh.count !== bonds.length * 2) {
      this.remove(this.bonds_mesh);
      this.bonds_mesh = this._get_bonds_mesh(bonds);
      this.add(this.bonds_mesh);
    }

    const matrix = new THREE.Matrix4();
    const dummy = new THREE.Object3D();
    const color = new THREE.Color();
    const distance = new THREE.Vector3();

    const particleA = new THREE.Object3D();
    const particleB = new THREE.Object3D();

    // create a flat list of all bonds using this.particles_mesh.count and getBondForParticle

    const allBonds = this._get_all_bonds(bonds);

    for (let instanceId = 0; instanceId < this.bonds_mesh.count; instanceId++) {
      const [particleAId, particleBId, bond_type] = allBonds[instanceId];
      this.particles_mesh.getMatrixAt(particleAId, matrix);
      matrix.decompose(particleA.position, particleA.rotation, particleA.scale);
      this.particles_mesh.getMatrixAt(particleBId, matrix);
      matrix.decompose(particleB.position, particleB.rotation, particleB.scale);

      this.bonds_mesh.getMatrixAt(instanceId, matrix);
      matrix.decompose(dummy.position, dummy.rotation, dummy.scale);

      dummy.position.copy(particleA.position);
      dummy.lookAt(particleB.position);

      distance.subVectors(particleA.position, particleB.position);
      const bondLength = distance.length();
      dummy.scale.set(1, 1, bondLength / 2);

      dummy.updateMatrix();
      this.bonds_mesh.setMatrixAt(instanceId, dummy.matrix);
      this.particles_mesh.getColorAt(particleAId, color);
      this.bonds_mesh.setColorAt(instanceId, color);
    }
    this.bonds_mesh.instanceMatrix.needsUpdate = true;
    this.bonds_mesh.instanceColor.needsUpdate = true;
  }

  step(frame) {
    if (frame !== undefined) {
      this.particle_cache = this.cache.get(frame);
    }
    if (this.particle_cache == null) {
      // nothing to display
    } else {
      this._updateParticles(this.particle_cache);

      if (this.show_bonds) {
        this._updateBonds(this.particle_cache.connectivity);
      }
    }
  }

  get_center(selection) {
    let selectedParticles = this.particle_cache;
    if (selection !== undefined && selection.length > 0) {
      selectedParticles = this.particle_cache.select(selection);
    }
    const center = new THREE.Vector3();
    center.copy(selectedParticles.getCenter());
    return center;
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
      if (document.activeElement === document.body) {
        if (event.repeat) {
          return;
        }
        if (event.isComposing || event.key === "i") {
          this.toggle();
        }
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
      // raycaster.camera = this.camera;
      const matrix = new THREE.Matrix4();
      const dummy = new THREE.Object3D();
      const direction = new THREE.Vector3();

      for (
        let instanceId = 0;
        instanceId < this.particlesGroup.particles_mesh.count;
        instanceId++
      ) {
        this.particlesGroup.particles_mesh.getMatrixAt(instanceId, matrix);
        matrix.decompose(dummy.position, dummy.rotation, dummy.scale);
        // combine all intersects from the center and top/bottom/left/right of the particle
        let visible = true;
        direction.copy(dummy.position).sub(this.camera.position).normalize();
        raycaster.set(this.camera.position, direction);
        const intersects = raycaster.intersectObjects(
          this.particlesGroup.children,
        );
        if (intersects.length > 0 && intersects[0].instanceId !== instanceId) {
          visible = false;
        }

        if (!visible) {
          // get the `${particle.name}-label`; object from this and remove it
          const label = this.getObjectByName(`${instanceId}-label`);
          this.remove(label);
        } else if (!this.getObjectByName(`${instanceId}-label`)) {
          // create a div with unicode f00d
          const text = document.createElement("div");
          const blank = "\u00A0";
          const line = "\u23AF";
          text.className = "label";
          if (this.label_offset === 0) {
            text.textContent = instanceId;
          } else if (this.label_offset > 0) {
            text.textContent = `${blank.repeat(
              this.label_offset * 2 + 6,
            )}\u00D7${line.repeat(this.label_offset)}${blank}${instanceId}`;
          } else {
            text.textContent = `${instanceId}${blank}${line.repeat(
              this.label_offset * -1,
            )}\u00D7${blank.repeat(this.label_offset * -2 + 6)}`;
          }
          text.style.fontSize = "20px";
          text.style.textShadow = "1px 1px 1px #000000";

          const label = new CSS2DObject(text);
          label.position.copy(dummy.position);
          label.name = `${instanceId}-label`;
          this.add(label);
        }
      }
    }
  }
}

export { ParticlesGroup, ParticleIndexGroup, CellGroup };
