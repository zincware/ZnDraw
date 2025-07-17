import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import * as znsocket from "znsocket";
import { client } from "../socket";
import { DEFAULT_ROOM_CONFIG } from "../types/room-config";
import { JMOL_COLORS, covalentRadii } from "./data";

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

export const useFrameConnection = (token: string) => {
    const framesCon = useMemo(() => {
        if (token === "") {
            return undefined;
        }
        // This single connection object will try the token-specific room first,
        // then automatically fall back to the default room if the key doesn't exist.
        return new znsocket.List({
            client,
            key: `room:${token}:frames`,
            fallback: "room:default:frames",
        });
    }, [token]);

    return { framesCon };
};

export const useFrameFetching = (token: string) => {
    const { framesCon } = useFrameConnection(token);

    const getFrameFromCon = useCallback(
        async (step: number) => {
            if (!framesCon) return null;
            return await framesCon.get(step);
        },
        [framesCon],
    );

    const getLengthFromCon = useCallback(async () => {
        if (!framesCon) return 0;
        return await framesCon.length();
    }, [framesCon]);

    return { getFrameFromCon, getLengthFromCon, framesCon };
};

export const setupFrames = (
    token: string,
    step: any,
    setCurrentFrame: any,
    currentFrame: any,
    setLength: any,
    setStep: any,
    frame_update: boolean,
    setIsFrameRendering?: (rendering: boolean) => void,
) => {
    const currentFrameUpdatedFromSocketRef = useRef(true);
    const [updateStepInPlace, setUpdateStepInPlace] = useState(0);
    const scaledRadii = useMemo(() => {
        const minRadius = Math.min(...covalentRadii);
        const maxRadius = Math.max(...covalentRadii);
        const range = maxRadius - minRadius;
        return covalentRadii.map((x: number) => (x - minRadius) / range + 0.3);
    }, [covalentRadii]);

    // Use the simplified fetching hook
    const { getFrameFromCon, getLengthFromCon, framesCon } =
        useFrameFetching(token);

    const setCurrentFrameFromObject = async (frame: any) => {
        const [positions, numbers, arrays, info, cell, constraints] =
            await Promise.all([
                frame.positions,
                frame.numbers,
                frame.arrays,
                frame.info,
                frame.cell,
                frame.constraints,
            ]);
        
        const [colors, radii, connectivity] = await Promise.all([
            arrays?.colors ?? Promise.resolve(null),
            arrays?.radii ?? Promise.resolve(null),
            info?.connectivity ?? Promise.resolve([]),
        ]);

        const resolvedPositions = positions.map(
            (position: [number, number, number]) =>
                new THREE.Vector3(position[0], position[1], position[2]),
        );

        const resolvedFrame = {
            ...frame,
            positions: resolvedPositions,
            numbers,
            connectivity: connectivity ?? [],
            cell: cell,
            constraints: constraints ?? [],
            arrays: {
                ...arrays,
                colors:
                    colors ??
                    numbers.map((x: number) => "#" + JMOL_COLORS[x].getHexString()),
                radii: radii ?? numbers.map((x: number) => scaledRadii[x]),
            },
        };

        setCurrentFrame(resolvedFrame);
        currentFrameUpdatedFromSocketRef.current = true;
        setIsFrameRendering?.(false);
    };

    useEffect(() => {
        if (framesCon === undefined) return;

        const handleRefresh = (x: any) => {
            if (x?.start > step && frame_update) {
                setStep(x.start);
            } else if (x?.start === step || x?.indices?.includes(step)) {
                setUpdateStepInPlace((prev) => prev + 1);
            }
        };

        framesCon.onRefresh(handleRefresh);

        return () => {
            framesCon.offRefresh?.(handleRefresh);
        };
    }, [framesCon, frame_update, step, setStep]);

    useEffect(() => {
        if (framesCon === undefined) {
            return;
        }

        const updateFrame = async () => {
            setIsFrameRendering?.(true);
            const length = await getLengthFromCon();
            setLength(length);
            if (step >= length) {
                setStep(Math.max(0, length - 1));
                setIsFrameRendering?.(false);
                return;
            }
            
            const frame = await getFrameFromCon(step || 0);

            if (frame === null) {
                console.warn("Frame ", step, " is null, retrying after 100ms...");
                setTimeout(updateFrame, 100);
                return;
            }
            setCurrentFrameFromObject(frame);
        };

        const debounceTimeout = setTimeout(updateFrame, 8);
        return () => clearTimeout(debounceTimeout);
    }, [step, updateStepInPlace, framesCon]);

    // Sending edits from ZnDraw to the server
    useEffect(() => {
        if (currentFrameUpdatedFromSocketRef.current === true) {
            currentFrameUpdatedFromSocketRef.current = false;
            return;
        }

        const updateCon = async () => {
            // Edits should only go to the token-specific room, not the fallback.
            // So we create a new connection without the fallback key for writing.
			// TODO: for edits, we need to ensure there exists a copy of the room
			// in the future we want to use segments!
            const con = new znsocket.List({
                client,
                key: `room:${token}:frames`,
            });
            if (!currentFrame.positions) {
                return;
            }
            
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

export const setupVectors = (
    token: string,
    step: any,
    setVectors: any,
    requestedProperties: string[] = [],
) => {
    // The hook now returns the simplified getFrameFromCon
    const { getFrameFromCon } = useFrameFetching(token);

    const setVectorsFromFrame = async (frame: any) => {
        if (!frame) {
            setVectors({});
            return;
        }

        try {
            const [calc, arrays] = await Promise.all([frame.calc, frame.arrays]);

            const [calcKeys, arraysKeys] = await Promise.all([
                calc.keys(),
                arrays.keys(),
            ]);

            const vectorProperties: { [key: string]: any } = {};
            const allKeys = new Set([...calcKeys, ...arraysKeys]);

            for (const property of requestedProperties) {
                if (!allKeys.has(property)) {
                    console.warn(`Requested property '${property}' not found in frame`);
                    continue;
                }

                let data;
                if (calcKeys.includes(property)) {
                    data = await calc[property];
                } else if (arraysKeys.includes(property)) {
                    data = await arrays[property];
                }

                if (
                    data &&
                    Array.isArray(data) &&
                    data.length > 0 &&
                    Array.isArray(data[0]) &&
                    data[0].length === 3
                ) {
                    vectorProperties[property] = data;
                } else {
                    console.warn(`Property '${property}' is not a valid 3D vector array`);
                }
            }
            setVectors(vectorProperties);
        } catch (error) {
            console.error("Error while resolving vectors:", error);
            setVectors({});
        }
    };

    useEffect(() => {
        const updateVectors = async () => {
            const frame = await getFrameFromCon(step || 0);
            if (frame === null) {
                console.warn(
                    "Frame ",
                    step,
                    " is null for vectors, retrying after 100ms...",
                );
                setTimeout(updateVectors, 100);
                return;
            }
            setVectorsFromFrame(frame);
        };

        const debounceTimeout = setTimeout(updateVectors, 8);
        return () => clearTimeout(debounceTimeout);
    }, [step, requestedProperties, getFrameFromCon]);
};
export const setupConfig = (token: string, setConfig: any) => {
	useEffect(() => {
		if (!token) return;

		const con = new znsocket.Dict({
			client: client,
			key: `room:${token}:config`,
		});

		// Get config sections from the DEFAULT_ROOM_CONFIG keys
		const configSections = Object.keys(DEFAULT_ROOM_CONFIG).filter(
			(key) => key !== "property",
		) as string[];
		const allConnections: any[] = [con];

		// Create connections for each config section
		const sectionConnections = configSections.map((section) => {
			const sectionCon = new znsocket.Dict({
				client: client,
				key: `room:${token}:config:${section}`,
			});
			allConnections.push(sectionCon);
			return { section, connection: sectionCon };
		});

		// Shared refresh handler that reloads the full config
		const refreshConfig = async (source: string) => {
			try {
				const result = await con.toObject();
				if (Object.keys(result).length > 0) {
					setConfig(result);
				} else {
					console.warn("Received empty config update");
				}
			} catch (error) {
				console.error(`Error refreshing config from ${source}:`, error);
			}
		};

		// Register refresh handlers for each section
		sectionConnections.forEach(({ section, connection }) => {
			connection.onRefresh(async () => {
				console.log(`${section} config section changed`);
				await refreshConfig(section);
			});
		});

		// Register main config refresh handler
		con.onRefresh(async () => {
			await refreshConfig("main");
		});

		// Initial load
		con
			.entries()
			.then((items: any) => con.toObject())
			.then((result) => {
				if (Object.keys(result).length > 0) {
					setConfig(result);
				}
			})
			.catch((error) => {
				console.error("Failed to load initial config:", error);
			});

		return () => {
			// Clean up all registered callbacks
			allConnections.forEach((connection) => {
				if (connection && typeof connection.offRefresh === "function") {
					connection.offRefresh();
				}
			});
		};
	}, [token]);
};
