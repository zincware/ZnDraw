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
      key: "room:" + token + ":points",
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
  const updateByRefreshRef = useRef(false);

  useEffect(() => {
    console.log("setting up selection");
    const con = new znsocket.Dict({
      client: client,
      key: "room:" + token + ":selection",
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
      console.log("updating selection");
      conInterfaceRef.current.set("grp-0", Array.from(selectedIds));
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
    con.get("grp-0").then((items: any) => {
      updateByRefreshRef.current = true;
      setStep(items);
    });

    con.onRefresh(async (x: any) => {
      const items = await con.get("grp-0");
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
      conInterfaceRef.current.set("grp-0", step);
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
      const items = Object.fromEntries(await con.entries());
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
      conInterface.set("position", cameraPosition.toArray());
      conInterface.set("target", orbitControlsTarget.toArray());
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
  currentFrame: any,
  setLength: any,
  setStep: any,
) => {
  const conInterfaceRef = useRef<typeof znsocket.List>(undefined);
  const defaultConInterfaceRef = useRef<typeof znsocket.List>(undefined);
  const useDefaultRoomRef = useRef(true);
  let [framesRequiringUpdate, setFramesRequiringUpdate] = useState(undefined);
  let currentFrameUpdatedFromSocketRef = useRef(true);

  const setCurrentFrameFromObject = (frame: any) => {
    frame = frame["value"];
    frame.positions = frame.positions.map(
      (position: [number, number, number]) =>
        new THREE.Vector3(position[0], position[1], position[2]),
    ) as THREE.Vector3[];
    setCurrentFrame(frame);
    currentFrameUpdatedFromSocketRef.current = true;
  };

  useEffect(() => {
    if (currentFrameUpdatedFromSocketRef.current === true) {
      currentFrameUpdatedFromSocketRef.current = false;
      return;
    }

    const updateCon = async () => {
      if (conInterfaceRef.current === undefined) {
        return;
      }
      if (!currentFrame.positions) {
        return;
      }
      // make the frame object serializable
      let newFrame = { ...currentFrame };
      newFrame.positions = newFrame.positions.map((x: THREE.Vector3) => [
        x.x,
        x.y,
        x.z,
      ]);
      // TODO: need to make a copy of the default room :/ Maybe when entering edit mode, send a message
      conInterfaceRef.current.set(parseInt(step) || 0, {
        value: newFrame,
        _type: "ase.Atoms",
      });
    };

    const debounceTimeout = setTimeout(updateCon, 500);

    return () => clearTimeout(debounceTimeout);
  }, [currentFrame]);

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
    // TODO with the edit mode, you need to check
    // if the room exists more often
    // because if can be updated from within
    // this instance and not trigger a refresh
    con.get(parseInt(step) || 0).then((frame: any) => {
      if (frame !== null) {
        useDefaultRoomRef.current = false;
        con.length().then((length: any) => {
          setLength(length);
        });
      } else {
        useDefaultRoomRef.current = true;
        defaultCon.length().then((length: any) => {
          setLength(length);
        });
      }
    });

    con.onRefresh(async (x: any) => {
      useDefaultRoomRef.current = false;
      con.length().then((length: any) => {
        setLength(length);
      });
      setFramesRequiringUpdate(x);
    });

    defaultCon.onRefresh(async (x: any) => {
      setFramesRequiringUpdate(x);

      defaultCon.length().then((length: any) => {
        setLength(length);
      });
    });

    conInterfaceRef.current = con;
    defaultConInterfaceRef.current = defaultCon;

    return () => {
      con.offRefresh();
      defaultCon.offRefresh();
    };
  }, [token]);

  useEffect(() => {
    if (
      conInterfaceRef.current === undefined &&
      defaultConInterfaceRef.current === undefined
    ) {
      return;
    }
    let currentInterface: znsocket.List = useDefaultRoomRef.current
      ? defaultConInterfaceRef.current
      : conInterfaceRef.current;
    if (framesRequiringUpdate !== undefined) {
      // cheap way out - we just update the current frame no matter what.
      currentInterface.length().then((length: any) => {
        setLength(length);
      });
      // the initial step is -1 so it queries the last frame
      currentInterface.get(parseInt(step) || 0).then((frame: any) => {
        if (frame === null) {
          currentInterface.length().then((length: any) => {
            setStep(length - 1);
          });
        } else {
          setCurrentFrameFromObject(frame);
        }
      });
      setFramesRequiringUpdate(undefined);
    }
  }, [step, framesRequiringUpdate]);

  useEffect(() => {
    if (
      conInterfaceRef.current === undefined &&
      defaultConInterfaceRef.current === undefined
    ) {
      return;
    }
    let currentInterface = useDefaultRoomRef.current
      ? defaultConInterfaceRef.current
      : conInterfaceRef.current;

    currentInterface.get(parseInt(step) || 0).then((frame: any) => {
      if (frame === null) {
        currentInterface.length().then((length: any) => {
          setStep(length - 1);
        });
      } else {
        setCurrentFrameFromObject(frame);
      }
    });
  }, [step]);
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
        geometries.push({
          ...value["value"]["data"],
          discriminator: value["value"]["class"],
        });
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
        geometries.push({
          ...value["value"]["data"],
          discriminator: value["value"]["class"],
        });
      }
      updateByRefreshRef.current = true;
      console.log(geometries);
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
      key: "room:" + token + ":config",
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
      if (Object.keys(result).length > 0) {
        setConfig(result);
      }
    });

    return () => {
      con.offRefresh();
    };
  }, [token]);
};
