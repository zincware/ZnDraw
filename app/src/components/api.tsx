import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import { JMOL_COLORS, covalentRadii } from "./data";
import * as znsocket from "znsocket";
import { client } from "../socket";

export function setupBookmarks(
	token: string,
	setBookmarks: any,
	bookmarks: any,
) {
	const conInterfaceRef = useRef(undefined);
	const updateByRefreshRef = useRef(false);
	useEffect(() => {
		if (token === "") {
			return;
		}
		const con = new znsocket.Dict({
			client: client,
			key: `room:${token}:bookmarks`,
		});

		// initial load
		con.entries().then((items: any) => {
			const result = Object.fromEntries(items);
			updateByRefreshRef.current = true;
			setBookmarks(result);
		});

		con.onRefresh(async (x: any) => {
			const items = await con.entries();
			const result = Object.fromEntries(items);
			console.log("bookmarks updated externally");
			updateByRefreshRef.current = true;
			setBookmarks(result);
		});
		conInterfaceRef.current = con;
		// setConInterface(con);

		return () => {
			con.offRefresh();
		};
	}, [token]);

	useEffect(() => {
		const updateCon = async () => {
			if (conInterfaceRef.current === undefined) {
				return;
			}
			console.log("updating bookmarks");
			conInterfaceRef.current.clear();
			await conInterfaceRef.current.update(bookmarks);
		};
		if (updateByRefreshRef.current) {
			updateByRefreshRef.current = false;
		} else {
			updateCon();
		}
	}, [bookmarks]);
}

export function setupPoints(token: string, setPoints: any, points: any) {
	const conInterfaceRef = useRef<any>(undefined);
	const updateByRefreshRef = useRef(true); // somewhere in the code the points are set initially
	useEffect(() => {
		const con = new znsocket.Dict({
			client: client,
			key: `room:${token}:points`,
		});

		// initial load
		con.get("grp-0").then((items: any) => {
			if (items === null) {
				return;
			}
			updateByRefreshRef.current = true;
			setPoints(items.map((x: number[]) => new THREE.Vector3(...x)));
		});

		con.onRefresh(async (x: any) => {
			const items = await con.get("grp-0");
			console.log("points updated externally");
			if (items === null) {
				return;
			}
			updateByRefreshRef.current = true;
			setPoints(items.map((x: number[]) => new THREE.Vector3(...x)));
		});

		conInterfaceRef.current = con;

		return () => {
			con.offRefresh();
		};
	}, [token]);

	useEffect(() => {
		const updateCon = async () => {
			if (conInterfaceRef.current === undefined) {
				return;
			}
			conInterfaceRef.current.clear();
			await conInterfaceRef.current.set(
				"grp-0",
				points.map((vec: THREE.Vector3) => [vec.x, vec.y, vec.z]),
			);
		};
		if (updateByRefreshRef.current) {
			updateByRefreshRef.current = false;
		} else {
			// Set a delay of 100ms before calling `updateCon`
			const debounceTimeout = setTimeout(updateCon, 100);
			// Clear the timeout if `points` changes within the 100ms
			return () => clearTimeout(debounceTimeout);
		}
	}, [points]);
}

export function setupSelection(
	token: string,
	setSelectedIds: (ids: Set<number>) => void,
	selectedIds: Set<number>,
) {
	const conInterfaceRef = useRef<any>(undefined);
	const updateByRefreshRef = useRef(true);

	useEffect(() => {
		if (token === "") {
			return;
		}
		const con = new znsocket.Dict({
			client: client,
			key: `room:${token}:selection`,
		});

		// Initial load
		con.get("grp-0").then((items: any) => {
			updateByRefreshRef.current = true;
			setSelectedIds(new Set(items));
		});

		// External updates handler
		con.onRefresh(async () => {
			const items = await con.get("grp-0");
			console.log("selection updated externally");
			updateByRefreshRef.current = true;
			setSelectedIds(new Set(items));
		});

		// Assign to ref instead of using useState
		conInterfaceRef.current = con;

		return () => {
			con.offRefresh();
		};
	}, [token]);

	useEffect(() => {
		const updateCon = async () => {
			if (!conInterfaceRef.current) return;

			// Update remote selection if change wasn't from external source
			conInterfaceRef.current.set("grp-0", Array.from(selectedIds));
		};

		if (updateByRefreshRef.current) {
			updateByRefreshRef.current = false;
		} else {
			const debounceTimeout = setTimeout(updateCon, 100);
			return () => clearTimeout(debounceTimeout);
		}
	}, [selectedIds]);
}

export function setupStep(
	token: string,
	setStep: (items: any) => void,
	step: any,
) {
	// Keep conInterfaceRef stable across renders
	const conInterfaceRef = useRef<any>(null);
	const updateByRefreshRef = useRef(false);

	// Memoize `con` to keep it consistent
	const con = useMemo(() => {
		return new znsocket.Dict({
			client: client,
			key: `room:${token}:step`,
		});
	}, [token]);

	useEffect(() => {
		// Load initial step only once
		const loadInitialStep = async () => {
			const items = await con.get("grp-0");
			console.log("initial step update");
			if (items !== null) {
				updateByRefreshRef.current = true;
				setStep(items);
			}
		};

		loadInitialStep();

		// Register external update listener
		con.onRefresh(async () => {
			const items = await con.get("grp-0");
			console.log("step updated externally");
			updateByRefreshRef.current = true;
			setStep(items);
		});

		conInterfaceRef.current = con;

		return () => {
			con.offRefresh();
		};
	}, [con, setStep]);

	const updateCon = useCallback(() => {
		if (conInterfaceRef.current) {
			conInterfaceRef.current.set("grp-0", step);
		}
	}, [step]);

	useEffect(() => {
		if (updateByRefreshRef.current) {
			// Reset flag if the step was updated externally
			updateByRefreshRef.current = false;
		} else {
			// Debounce the update to avoid rapid consecutive updates
			const debounceTimeout = setTimeout(updateCon, 100);
			return () => clearTimeout(debounceTimeout);
		}
	}, [step, updateCon]);
}
export function setupCamera(
	token: string,
	cameraAndControls: { camera: any; target: any },
	setCameraAndControls: (cameraAndControls: {
		camera: any;
		target: any;
	}) => void,
	synchronize_camera: boolean,
) {
	// let [conInterface, setConInterface]: any = useState(undefined);

	const conInterface = useMemo(() => {
		return new znsocket.Dict({
			client: client,
			key: `room:${token}:camera`,
		});
	}, [token]);

	const updateByRefreshRef = useRef(true);

	useEffect(() => {
		if (synchronize_camera) {
			conInterface.onRefresh(async (x: any) => {
				const items = Object.fromEntries(await conInterface.entries());
				console.log("camera updated externally");
				setCameraAndControls({
					camera: new THREE.Vector3(...items.position),
					target: new THREE.Vector3(...items.target),
				});

				// TODO
				updateByRefreshRef.current = true;
			});
		}

		return () => {
			conInterface.offRefresh();
		};
	}, [token, synchronize_camera, conInterface]);

	// initial load of the camera and target
	useEffect(() => {
		const loadInitialCamera = async () => {
			const items = Object.fromEntries(await conInterface.entries());
			if (
				items === null ||
				items.position === undefined ||
				items.target === undefined
			) {
				return;
			}
			setCameraAndControls({
				camera: new THREE.Vector3(...items.position),
				target: new THREE.Vector3(...items.target),
			});
		};

		loadInitialCamera();
	}, [conInterface, setCameraAndControls]);

	useEffect(() => {
		const updateCon = async () => {
			if (updateByRefreshRef.current) {
				// this is triggered twice because of camera and orbit target
				// TODO: this is not reliable yet!
				updateByRefreshRef.current = false;
				return;
			}
			conInterface.update({
				position: cameraAndControls.camera.toArray(),
				target: cameraAndControls.target.toArray(),
			});
		};

		// Set a delay of 100ms before calling `updateCon`
		const debounceTimeout = setTimeout(updateCon, 100);

		// Clear the timeout if `points` changes within the 100ms
		return () => clearTimeout(debounceTimeout);
	}, [cameraAndControls, conInterface]);
}

export const setupFrames = (
	token: string,
	step: any,
	setCurrentFrame: any,
	currentFrame: any,
	setLength: any,
	setStep: any,
	frame_update: boolean,
) => {
	const currentFrameUpdatedFromSocketRef = useRef(true);
	const customRoomAvailRef = useRef(false); // Track whether listening to the default room
	const [updateStepInPlace, setUpdateStepInPlace] = useState(0);
	const scaledRadii = useMemo(() => {
		const minRadius = Math.min(...covalentRadii);
		const maxRadius = Math.max(...covalentRadii);
		const range = maxRadius - minRadius;
		return covalentRadii.map((x: number) => (x - minRadius) / range + 0.3);
	}, [covalentRadii]);

	const setCurrentFrameFromObject = async (frame: any) => {
		// Await top-level fields
		const [positions, numbers, arrays, connectivity, cell] = await Promise.all([
			frame.positions,
			frame.numbers,
			frame.arrays,
			frame.connectivity,
			frame.cell,
		]);

		// TODO: vectors are missing!
		// request them seperately, only if needed

		// Await nested arrays.colors and arrays.radii
		const [colors, radii] = await Promise.all([
			arrays?.colors ?? Promise.resolve(null),
			arrays?.radii ?? Promise.resolve(null),
		]);

		// Convert positions to THREE.Vector3
		const resolvedPositions = positions.map(
			(position: [number, number, number]) =>
				new THREE.Vector3(position[0], position[1], position[2]),
		);

		// Build full resolved frame
		const resolvedFrame = {
			...frame,
			positions: resolvedPositions,
			numbers,
			connectivity: connectivity ?? [],
			cell: cell,
			arrays: {
				...arrays,
				colors: colors ??
					numbers.map(
						(x: number) => "#" + JMOL_COLORS[x].getHexString(),
					),
				radii: radii ?? numbers.map((x: number) => scaledRadii[x]),
			},
		};

		console.log("Frame resolved from socket:", resolvedFrame);

		setCurrentFrame(resolvedFrame);
		currentFrameUpdatedFromSocketRef.current = true;
	};

	const currentRoomCon = useMemo(() => {
		if (token === "") {
			return undefined;
		}
		return new znsocket.List({
			client,
			key: `room:${token}:frames`,
		});
	}, [token]);

	const defaultRoomCon = useMemo(() => {
		return new znsocket.List({
			client,
			key: "room:default:frames",
		});
	}, []);

	// setup onRefresh and initial load
	// default room
	useEffect(() => {
		if (defaultRoomCon === undefined) return;

		defaultRoomCon.onRefresh(async (x: any) => {
			if (x?.start > step && frame_update) {
				setStep(x.start);
			} else if (x?.start === step) {
				setUpdateStepInPlace((prev) => prev + 1);
			} else if (x?.indices?.includes(step)) {
				setUpdateStepInPlace((prev) => prev + 1);
			}
		});

		return () => {
			defaultRoomCon.offRefresh?.();
		};
	}, [defaultRoomCon, frame_update, step]);

	// custom room
	useEffect(() => {
		if (currentRoomCon === undefined) return;

		currentRoomCon.onRefresh(async (x: any) => {
			if (x?.start > step && frame_update) {
				setStep(x.start);
			} else if (x?.start === step) {
				setUpdateStepInPlace((prev) => prev + 1);
			} else if (x?.indices?.includes(step)) {
				setUpdateStepInPlace((prev) => prev + 1);
			}
		});

		return () => {
			currentRoomCon.offRefresh?.();
		};
	}, [currentRoomCon, frame_update, step]);

	const getFrameFromCon = useCallback(
		async (step: number) => {
			let fetchedFrame = await currentRoomCon.get(step);
			if (!customRoomAvailRef.current && fetchedFrame === null) {
				return await defaultRoomCon.get(step);
			}
			customRoomAvailRef.current = true;
			return fetchedFrame;
		},
		[currentRoomCon, defaultRoomCon],
	);

	const getLengthFromCon = useCallback(async () => {
		let length = await currentRoomCon.length();
		if (!customRoomAvailRef.current && length === 0) {
			return await defaultRoomCon.length();
		}
		customRoomAvailRef.current = true;
		return length;
	}, [currentRoomCon, defaultRoomCon]);

	useEffect(() => {
		if (defaultRoomCon === undefined || currentRoomCon === undefined) {
			return;
		}

		const updateFrame = async () => {
			// first we check the length
			const length = await getLengthFromCon();
			setLength(length);
			if (step >= length) {
				setStep(Math.max(0, length - 1));
				return;
			}
			// now we request the frame

			const frame = await getFrameFromCon(step || 0);
			

			if (frame === null) {
				console.warn("Frame ", step, " is null, retrying after 100ms...");
				setTimeout(updateFrame, 100); // Retry after 100 ms
				return;
			}
			setCurrentFrameFromObject(frame); // Process the valid frame
		};

		// debounce by 8ms (roughly 120fps)
		// this will smooth scrubbing through the frames
		// and avoid unnecessary updates
		const debounceTimeout = setTimeout(updateFrame, 8);
		return () => {
			clearTimeout(debounceTimeout);
		};
	}, [step, updateStepInPlace, defaultRoomCon, currentRoomCon]);

	// Sending edits from ZnDraw to the server
	useEffect(() => {
		if (currentFrameUpdatedFromSocketRef.current === true) {
			currentFrameUpdatedFromSocketRef.current = false;
			return;
		}

		const updateCon = async () => {
			const con = new znsocket.List({
				client,
				key: `room:${token}:frames`,
			});
			if (!currentFrame.positions) {
				return;
			}
			// Make the frame object serializable
			const newFrame = { ...currentFrame };
			newFrame.positions = newFrame.positions.map((x: THREE.Vector3) => [
				x.x,
				x.y,
				x.z,
			]);
			con.set(Number.parseInt(step) || 0, {
				value: newFrame,
				_type: "ase.Atoms",
			});
		};

		// Use a debounce mechanism for updates
		const debounceTimeout = setTimeout(updateCon, 100);

		return () => clearTimeout(debounceTimeout);
	}, [currentFrame]);
};

export const setupGeometries = (
	token: string,
	setGeometries: any,
	geometries: any,
) => {
	const conInterfaceRef = useRef<typeof znsocket.List>(undefined);
	const updateByRefreshRef = useRef<boolean>(false);

	useEffect(() => {
		const con = new znsocket.List({
			client: client,
			key: `room:${token}:geometries`,
		});

		// initial load
		const loadGeometries = async () => {
			const geometries = [];
			for await (const value of con) {
				geometries.push({
					...value.value.data,
					discriminator: value.value.class,
				});
			}
			updateByRefreshRef.current = true;
			setGeometries(geometries);
		};
		loadGeometries();

		con.onRefresh(async (x: any) => {
			// console.log(x); one could only update the indices or simply all. We go for all
			// async iterate the conInterface
			const geometries = [];
			console.log("geometries updated externally");
			for await (const value of con) {
				// console.log(value["value"]["data"]); // Process each value as it becomes available
				geometries.push({
					...value.value.data,
					discriminator: value.value.class,
				});
			}
			updateByRefreshRef.current = true;
			setGeometries(geometries);
		});

		conInterfaceRef.current = con;

		return () => {
			con.offRefresh();
		};
	}, [token]);

	useEffect(() => {
		const updateCon = async () => {
			if (conInterfaceRef.current === undefined) {
				return;
			}
			conInterfaceRef.current.clear();
			for (let i = 0; i < geometries.length; i++) {
				await conInterfaceRef.current.push({ data: geometries[i] });
			}
		};

		if (updateByRefreshRef.current) {
			updateByRefreshRef.current = false;
		} else {
			updateCon();
		}
	}, [geometries]);
};

export const setupFigures = (token: string, setUpdatedPlotsList: any) => {
	useEffect(() => {
		const con = new znsocket.Dict({
			client: client,
			key: `room:${token}:figures`,
		});

		con.onRefresh(async (x: any) => {
			setUpdatedPlotsList(x.keys);
		});
		return () => {
			con.offRefresh();
		};
	}, [token]);
};

export const setupMessages = (
	token: string,
	setMessages: any,
	messages: any[],
) => {
	const conInterfaceRef = useRef<typeof znsocket.List>(undefined);
	const updateByRefreshRef = useRef<boolean>(false);

	useEffect(() => {
		const con = new znsocket.List({
			client: client,
			key: `room:${token}:chat`,
		});

		// initial load
		const loadMessages = async () => {
			const messages = [];
			for await (const value of con) {
				messages.push(value);
			}
			updateByRefreshRef.current = true;
			setMessages(messages);
		};
		loadMessages();

		con.onRefresh(async (x: any) => {
			const messages = [];
			for await (const value of con) {
				messages.push(value);
			}
			updateByRefreshRef.current = true;
			setMessages(messages);
		});

		conInterfaceRef.current = con;

		return () => {
			con.offRefresh();
		};
	}, [token]);

	useEffect(() => {
		const updateCon = async () => {
			if (conInterfaceRef.current === undefined) {
				return;
			}
			conInterfaceRef.current.clear();
			for (let i = 0; i < messages.length; i++) {
				await conInterfaceRef.current.push(messages[i]);
			}
		};

		if (updateByRefreshRef.current) {
			updateByRefreshRef.current = false;
		} else {
			updateCon();
		}
	}, [messages]);
};

export const setupConfig = (token: string, setConfig: any) => {
	useEffect(() => {
		const con = new znsocket.Dict({
			client: client,
			key: `room:${token}:config`,
		});

		// initial load
		con
			.entries()
			.then((items: any) => con.toObject())
			.then((result) => {
				if (Object.keys(result).length > 0) {
					setConfig(result);
				}
			})
			.catch((error) => {
				console.error("Failed to load config:", error);
			});

		con.onRefresh(async (x: any) => {
			const result = await con.toObject();
			console.log("config updated externally");
			if (Object.keys(result).length > 0) {
				setConfig(result);
			}
		});

		return () => {
			con.offRefresh();
		};
	}, [token]);
};
