import * as THREE from 'three';

export const atomsGroup = new THREE.Group();

export const bondsGroup = new THREE.Group();
export const bondsGroup_1 = new THREE.Group();
export const bondsGroup_2 = new THREE.Group();

export const materials = {
	"MeshBasicMaterial": new THREE.MeshBasicMaterial({ color: "#ffa500" }),
	"MeshLambertMaterial": new THREE.MeshLambertMaterial({ color: "#ffa500" }),
	"MeshMatcapMaterial": new THREE.MeshMatcapMaterial({ color: "#ffa500" }),
	"MeshPhongMaterial": new THREE.MeshPhongMaterial({ color: "#ffa500" }),
	"MeshPhysicalMaterial": new THREE.MeshPhysicalMaterial({ color: "#ffa500" }),
	"MeshStandardMaterial": new THREE.MeshStandardMaterial({ color: "#ffa500" }),
	"MeshToonMaterial": new THREE.MeshToonMaterial({ color: "#ffa500" }),

};

let node1 = new THREE.Vector3();
let node2 = new THREE.Vector3();

// TODO do not create new variables every time but reuse!
function halfCylinderGeometry(pointX, pointY, config) {
	// Make the geometry (of "direction" length)
	var direction = new THREE.Vector3().subVectors(pointY, pointX);

	var geometry = new THREE.CylinderGeometry(0.15 * config["bond_size"], 0.15 * config["bond_size"], direction.length() / 2, config["resolution"] * 2);
	// // shift it so one end rests on the origin
	geometry.applyMatrix4(new THREE.Matrix4().makeTranslation(0, direction.length() / 4, 0));
	// // rotate it the right way for lookAt to work
	geometry.applyMatrix4(new THREE.Matrix4().makeRotationX(THREE.MathUtils.degToRad(90)));

	return geometry;

}

function halfCylinderMesh(pointX, pointY, material, config) {
	// // Make a mesh with the geometry
	var geometry = halfCylinderGeometry(pointX, pointY, config);
	var mesh = new THREE.Mesh(geometry, material);
	// alignBetweenVectors(pointX, pointY, mesh);
	// // Position it where we want
	mesh.position.copy(pointX);
	// // And make it point to where we want
	mesh.lookAt(pointY);

	return mesh;
}

function addAtom(item, config) {
	const geometry = new THREE.SphereGeometry(item["radius"] * config["sphere_size"], config["resolution"] * 4, config["resolution"] * 2);

	const particle = new THREE.Mesh(geometry, materials[document.getElementById('materialSelect').value].clone());
	atomsGroup.add(particle);
	particle.material.color.set(item["color"]);
	particle.material.wireframe = document.getElementById('wireframe').checked;

	particle.position.set(...item["position"]);
	particle.userData["id"] = item["id"];
	particle.userData["color"] = item["color"];
	particle.userData["bond_ids"] = [];
}

function addBond(item, config) {
	atomsGroup.children[item[0]].getWorldPosition(node1);
	atomsGroup.children[item[1]].getWorldPosition(node2);

	const bond_1 = halfCylinderMesh(node1, node2, atomsGroup.children[item[0]].material, config);
	const bond_2 = halfCylinderMesh(node2, node1, atomsGroup.children[item[1]].material, config);

	bond_1.userData["atom_id"] = item[0];
	bond_2.userData["atom_id"] = item[1];

	atomsGroup.children[item[0]].userData["bond_ids"].push(bond_1.id);
	atomsGroup.children[item[0]].userData["bond_ids"].push(bond_2.id);
	atomsGroup.children[item[1]].userData["bond_ids"].push(bond_1.id);
	atomsGroup.children[item[1]].userData["bond_ids"].push(bond_2.id);


	bondsGroup_1.add(bond_1);
	bondsGroup_2.add(bond_2);

}

// TODO rename clean AtomsBonds
// this does not remove the scene!
export function cleanScene(scene) {
	while (atomsGroup.children.length > 0) {
		scene.remove(atomsGroup.children.shift());
	};
	while (bondsGroup_1.children.length > 0) {
		scene.remove(bondsGroup_1.children.shift());
	};
	while (bondsGroup_2.children.length > 0) {
		scene.remove(bondsGroup_2.children.shift());
	};
}

export function drawAtoms(atoms, bonds, config, scene) {
	cleanScene(scene);

	atoms.forEach(function (item, index) {
		// console.log("Adding item " + index + " to scene(" + item + ")");
		addAtom(item, config);
	});

	scene.add(atomsGroup);
	if (config["bond_size"] > 0) {

		bonds.forEach(function (item, index) {
			// console.log("Adding item " + index + " to scene(" + item + ")");
			addBond(item, config);

		});

		bondsGroup.add(bondsGroup_1);
		bondsGroup.add(bondsGroup_2);
		scene.add(bondsGroup);

	}

}