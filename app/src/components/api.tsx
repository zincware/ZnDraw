import { useEffect, useState } from "react";
import { client } from "../socket";
import * as znsocket from "znsocket";
import * as THREE from "three";

export function setupBookmarks(
  token: string,
  setBookmarks: any,
  bookmarks: any,
) {
  let [conInterface, setConInterface]: any = useState(undefined);
  useEffect(() => {
    const con = new znsocket.Dict({
      client: client,
      key: "room:" + token + ":bookmarks",
    });

    con.onRefresh(async (x: any) => {
      const items = await con.items();
      const result = Object.fromEntries(items);
      console.log("bookmarks updated externally");
      setBookmarks(result);
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
      // TODO use bookmarksInterface.update
      conInterface.clear();
      for (let key in bookmarks) {
        await conInterface.setitem(parseInt(key), bookmarks[key]);
      }
    };
    updateCon();
  }, [bookmarks, conInterface]);
}

export function setupPoints(token: string, setPoints: any, points: any) {
  let [conInterface, setConInterface]: any = useState(undefined);
  useEffect(() => {
    const con = new znsocket.Dict({
      client: client,
      key: "room:" + token + ":points",
    });

    con.onRefresh(async (x: any) => {
      const items = await con.getitem("0");
      const result = Object.fromEntries(items);
      setPoints(result["0"].map((x: number[]) => new THREE.Vector3(...x)));
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
      conInterface.clear();
      await conInterface.setitem(
        0,
        points.map((vec: THREE.Vector3) => [vec.x, vec.y, vec.z]),
      );
      // TODO: set segments
    };

    // Set a delay of 100ms before calling `updateCon`
    const debounceTimeout = setTimeout(updateCon, 100);

    // Clear the timeout if `points` changes within the 100ms
    return () => clearTimeout(debounceTimeout);
  }, [points, conInterface]);
}

export function setupSelection(
  token: string,
  setSelection: any,
  selection: any,
) {
  let [conInterface, setConInterface]: any = useState(undefined);
  useEffect(() => {
    const con = new znsocket.Dict({
      client: client,
      key: "room:" + token + ":selection",
    });

    con.onRefresh(async (x: any) => {
      const items = await con.getitem(0);
      console.log("selection updated externally");
      setSelection(new Set(items));
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
      conInterface.setitem(0, Array.from(selection));
    };
    updateCon();
  }, [selection, conInterface]);
}

export function setupStep(token: string, setStep: any, step: any) {
  let [conInterface, setConInterface]: any = useState(undefined);
  useEffect(() => {
    const con = new znsocket.Dict({
      client: client,
      key: "room:" + token + ":step",
    });

    con.onRefresh(async (x: any) => {
      const items = await con.getitem(0);
      console.log("step updated externally");
      setStep(items);
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
      conInterface.setitem(0, step);
    };
    // Set a delay of 100ms before calling `updateCon`
    const debounceTimeout = setTimeout(updateCon, 100);

    // Clear the timeout if `points` changes within the 100ms
    return () => clearTimeout(debounceTimeout);
  }, [step, conInterface]);
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
  useEffect(() => {
    const con = new znsocket.Dict({
      client: client,
      key: "room:" + token + ":camera",
    });

    con.onRefresh(async (x: any) => {
      const items = Object.fromEntries(await con.items());
      console.log(items);
      setCameraPosition(new THREE.Vector3(...items.position));
      setOrbitControlsTarget(new THREE.Vector3(...items.target));

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

export const setupFrames = (token: string, step: any, setCurrentFrame: any, setLength: any) => {
  let [conInterface, setConInterface]: any = useState(undefined);
  let [useDefaultRoom, setUseDefaultRoom] = useState(true);

  useEffect(() => {
    const con = new znsocket.List({
      client: client,
      key: "room:" + token + ":frames",
    });
    // initially check if the room exists
    con.getitem(parseInt(step)).then((frame: any) => {
      if (frame === null) {
        setUseDefaultRoom(true);
      }
    });

    con.onRefresh(async (x: any) => {
      setUseDefaultRoom(false);

      con.len().then((length: any) => {setLength(length)});

      if (x["start"] && x["stop"]) {
        if (parseInt(step) > x["start"] || parseInt(step) < x["stop"]) {
          // update in place, also stop might not be defined
        }
      } else if (x["indices"]) {
        if (x["indices"].includes(parseInt(step))) {
          // update in place
        }
      }
    });

    setConInterface(con);
    if (useDefaultRoom){
      const defaultCon = new znsocket.List({
        client: client,
        key: "room:default:frames",
      });
      defaultCon.len().then((length: any) => {setLength(length)});
    } else {
      con.len().then((length: any) => {setLength(length)});
    }

    return () => {
      con.offRefresh();
    };
  }, [token]);

  useEffect(() => {
    if (conInterface === undefined) {
      return;
    }
    let currentInterface = conInterface;
    if (useDefaultRoom) {
      currentInterface = new znsocket.List({
        client: client,
        key: "room:default:frames",
      });
    }

    currentInterface.getitem(parseInt(step)).then((frame: any) => {
      frame = frame["value"];
      frame.positions = frame.positions.map(
        (position: [number, number, number]) =>
          new THREE.Vector3(position[0], position[1], position[2]),
      ) as THREE.Vector3[];
      setCurrentFrame(frame);
    });

  }, [conInterface, step, useDefaultRoom]);
};