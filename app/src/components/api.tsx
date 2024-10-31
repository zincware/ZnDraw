import { useEffect, useState, useRef } from "react";
import { client } from "../socket";
import * as znsocket from "znsocket";
import * as THREE from "three";

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
      key: "room:" + token + ":bookmarks",
    });

    // initial load
    con.items().then((items: any) => {
      const result = Object.fromEntries(items);
      updateByRefreshRef.current = true;
      setBookmarks(result);
    });

    con.onRefresh(async (x: any) => {
      const items = await con.items();
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
      key: "room:" + token + ":points",
    });

    // initial load
    con.getitem(0).then((items: any) => {
      if (items === null) {
        return;
      }
      updateByRefreshRef.current = true;
      setPoints(items.map((x: number[]) => new THREE.Vector3(...x)));
    });

    con.onRefresh(async (x: any) => {
      const items = await con.getitem(0);
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
      await conInterfaceRef.current.setitem(
        0,
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
  const updateByRefreshRef = useRef(false);

  useEffect(() => {
    console.log("setting up selection");
    const con = new znsocket.Dict({
      client: client,
      key: "room:" + token + ":selection",
    });

    // Initial load
    con.getitem(0).then((items: any) => {
      updateByRefreshRef.current = true;
      setSelectedIds(new Set(items));
    });

    // External updates handler
    con.onRefresh(async () => {
      const items = await con.getitem(0);
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
      console.log("updating selection");
      conInterfaceRef.current.setitem(0, Array.from(selectedIds));
    };

    if (updateByRefreshRef.current) {
      updateByRefreshRef.current = false;
    } else {
      updateCon();
    }
  }, [selectedIds]);
}

export function setupStep(token: string, setStep: any, step: any) {
  const conInterfaceRef = useRef<any>(undefined);
  const updateByRefreshRef = useRef(false);

  useEffect(() => {
    const con = new znsocket.Dict({
      client: client,
      key: "room:" + token + ":step",
    });

    // Initial load
    con.getitem(0).then((items: any) => {
      updateByRefreshRef.current = true;
      setStep(items);
    });

    con.onRefresh(async (x: any) => {
      const items = await con.getitem(0);
      console.log("step updated externally");
      updateByRefreshRef.current = true;
      setStep(items);
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
      conInterfaceRef.current.setitem(0, step);
    };

    if (updateByRefreshRef.current) {
      updateByRefreshRef.current = false;
    } else {
      // Set a delay of 100ms before calling `updateCon`
      const debounceTimeout = setTimeout(updateCon, 100);

      // Clear the timeout if `points` changes within the 100ms
      return () => clearTimeout(debounceTimeout);
    }
  }, [step]);
}

export function setupCamera(
  token: string,
  cameraPosition: any,
  orbitControlsTarget: any,
  setCameraPosition: any,
  setOrbitControlsTarget: any,
  controlsRef: any,
  cameraRef: any,
) {
  let [conInterface, setConInterface]: any = useState(undefined);
  const updateByRefreshRef = useRef(0);
  useEffect(() => {
    const con = new znsocket.Dict({
      client: client,
      key: "room:" + token + ":camera",
    });

    con.onRefresh(async (x: any) => {
      const items = Object.fromEntries(await con.items());
      console.log("camera updated externally");
      setCameraPosition(new THREE.Vector3(...items.position));
      setOrbitControlsTarget(new THREE.Vector3(...items.target));

      updateByRefreshRef.current = 2;

      if (controlsRef.current && cameraRef.current) {
        controlsRef.current.enabled = false;
        cameraRef.current.position.set(...items.position);
        controlsRef.current.update();
        controlsRef.current.enabled = true;
      }
    });

    setConInterface(con);

    return () => {
      con.offRefresh();
    };
  }, [token]);

  useEffect(() => {
    const updateCon = async () => {
      if (conInterface === undefined) {
        return;
      }
      if (updateByRefreshRef.current > 0) {
        // this is triggered twice because of camera and orbit target
        // TODO: this is not reliable yet!
        updateByRefreshRef.current -= 1;
        return;
      }
      // TODO: use conInterface.update
      conInterface.setitem("position", cameraPosition.toArray());
      conInterface.setitem("target", orbitControlsTarget.toArray());
    };

    // Set a delay of 100ms before calling `updateCon`
    const debounceTimeout = setTimeout(updateCon, 100);

    // Clear the timeout if `points` changes within the 100ms
    return () => clearTimeout(debounceTimeout);
  }, [cameraPosition, orbitControlsTarget, conInterface]);
}

export const setupFrames = (
  token: string,
  step: any,
  setCurrentFrame: any,
  setLength: any,
  setStep: any,
) => {
  let [conInterface, setConInterface]: any = useState(undefined); // TODO: useRef instead
  let [defaultConInterface, setDefaultConInterface]: any = useState(undefined); // TODO: useRef instead
  let [useDefaultRoom, setUseDefaultRoom] = useState(true); // TODO: useRef instead?
  let [framesRequiringUpdate, setFramesRequiringUpdate] = useState(undefined);

  const setCurrentFrameFromObject = (frame: any) => {
    console.log("setting current frame");
    frame = frame["value"];
    frame.positions = frame.positions.map(
      (position: [number, number, number]) =>
        new THREE.Vector3(position[0], position[1], position[2]),
    ) as THREE.Vector3[];
    setCurrentFrame(frame);
  };

  useEffect(() => {
    const con = new znsocket.List({
      client: client,
      key: "room:" + token + ":frames",
    });

    const defaultCon = new znsocket.List({
      client: client,
      key: "room:default:frames",
    });

    // initially check if the room exists
    con.getitem(parseInt(step) || 0).then((frame: any) => {
      if (frame !== null) {
        setUseDefaultRoom(false);
        con.len().then((length: any) => {
          setLength(length);
        });
      } else {
        setUseDefaultRoom(true);
        defaultCon.len().then((length: any) => {
          setLength(length);
        });
      }
    });

    con.onRefresh(async (x: any) => {
      setUseDefaultRoom(false);
      con.len().then((length: any) => {
        setLength(length);
      });
      setFramesRequiringUpdate(x);
    });

    defaultCon.onRefresh(async (x: any) => {
      setFramesRequiringUpdate(x);

      defaultCon.len().then((length: any) => {
        setLength(length);
      });
    });

    setConInterface(con);
    setDefaultConInterface(defaultCon);

    return () => {
      con.offRefresh();
      defaultCon.offRefresh();
    };
  }, [token]);

  useEffect(() => {
    if (conInterface === undefined && defaultConInterface === undefined) {
      return;
    }
    let currentInterface = useDefaultRoom ? defaultConInterface : conInterface;
    if (framesRequiringUpdate !== undefined) {
      // cheap way out - we just update the current frame no matter what.
      currentInterface.len().then((length: any) => {
        setLength(length);
      });
      currentInterface.getitem(parseInt(step) || 0).then((frame: any) => {
        console.log(frame);
        if (frame === null) {
          currentInterface.len().then((length: any) => {
            setStep(length - 1);
          });
        } else {
          setCurrentFrameFromObject(frame);
        }
      });
      setFramesRequiringUpdate(undefined);
    }
  }, [
    conInterface,
    step,
    useDefaultRoom,
    defaultConInterface,
    framesRequiringUpdate,
  ]);

  useEffect(() => {
    if (conInterface === undefined && defaultConInterface === undefined) {
      return;
    }
    let currentInterface = useDefaultRoom ? defaultConInterface : conInterface;

    currentInterface.getitem(parseInt(step) || 0).then((frame: any) => {
      if (frame === null) {
        currentInterface.len().then((length: any) => {
          setStep(length - 1);
        });
      } else {
        setCurrentFrameFromObject(frame);
      }
    });
  }, [conInterface, step, useDefaultRoom, defaultConInterface]);
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
      key: "room:" + token + ":geometries",
    });

    // initial load
    const loadGeometries = async () => {
      let geometries = [];
      for await (const value of con) {
        geometries.push(value["value"]["data"]);
      }
      updateByRefreshRef.current = true;
      setGeometries(geometries);
    };
    loadGeometries();

    con.onRefresh(async (x: any) => {
      // console.log(x); one could only update the indices or simply all. We go for all
      // async iterate the conInterface
      let geometries = [];
      console.log("geometries updated externally");
      for await (const value of con) {
        // console.log(value["value"]["data"]); // Process each value as it becomes available
        geometries.push(value["value"]["data"]);
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
        await conInterfaceRef.current.append({ data: geometries[i] });
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
      key: "room:" + token + ":figures",
    });

    con.onRefresh(async (x: any) => {
      setUpdatedPlotsList(x["keys"]);
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
      key: "room:" + token + ":chat",
    });

    // initial load
    const loadMessages = async () => {
      let messages = [];
      for await (const value of con) {
        messages.push(value);
      }
      updateByRefreshRef.current = true;
      setMessages(messages);
    };
    loadMessages();

    con.onRefresh(async (x: any) => {
      let messages = [];
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
        await conInterfaceRef.current.append(messages[i]);
      }
    };

    if (updateByRefreshRef.current) {
      updateByRefreshRef.current = false;
    } else {
      updateCon();
    }
  }, [messages]);
};
