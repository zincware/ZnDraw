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
            } else {
                const particle = new THREE.Mesh(
                    new THREE.SphereGeometry(item.radius, 32, 32),
                    new THREE.MeshBasicMaterial({ color: item.color })
                );
                particle.position.set(item.x, item.y, item.z);
                particle.name = item.id;
                particleGroup.add(particle);
         }
         });
    }

    return particleGroup;
}