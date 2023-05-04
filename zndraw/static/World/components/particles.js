import * as THREE from 'three';


export function createParticleGroup() {
    const particleGroup = new THREE.Group();

    particleGroup.tick = (data) => {
        if (data == null) {
            return;
        }
        const particles = data["particles"];
        const bonds = data["bonds"];
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

        // Update bonds
        bonds.forEach((item) => {
            const geometry = new THREE.BufferGeometry().setFromPoints(
                [
                    particleGroup.getObjectByName(item[0]).position,
                    particleGroup.getObjectByName(item[1]).position,
                ]
            );
            const material = new THREE.LineBasicMaterial({
                color: 0xff00ff
            });
            const line = new THREE.Line( geometry, material );
            particleGroup.add(line);
            // console.log(item);
        });

    }

    return particleGroup;
}