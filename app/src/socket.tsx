import { useEffect, useRef } from "react";
import { Manager } from "socket.io-client";
import { createClient } from "znsocket";

import * as THREE from "three";

function setupIO() {
	const basePath = import.meta.env.BASE_URL || "/";
	let manager;

	if (basePath === "/") {
		// manager = new Manager("http://localhost:1235"); // for local development
		manager = new Manager(window.location.origin); // for production
	} else {
		manager = new Manager(window.location.origin, {
			path: `${basePath}socket.io`,
		});
	}
	console.log("manager", manager);
	return {
		socket: manager.socket("/"),
		client: createClient({ socket: manager.socket("/znsocket") }),
	};
}
export const { socket, client } = setupIO();
