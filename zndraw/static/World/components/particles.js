import * as THREE from 'three';


export function createParticleGroup() {
    const particleGroup = new THREE.Group();
    
    particleGroup.tick = (data) => {
        if (data == null) {
            return;
        }
        data.forEach((item) => { 
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
                if (data.filter((item) => item.id === particle.name).length === 0) {
                    particleGroup.remove(particle);
                }
            }
        );

    }

    return particleGroup;
}