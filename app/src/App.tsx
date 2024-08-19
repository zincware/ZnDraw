import React, { useState, useEffect, useRef } from "react";
import {
  socket,
  sendStep,
  sendSelection,
  sendBookmarks,
  sendCamera,
  sendPoints,
} from "./socket";
import HeadBar from "./components/headbar";
import Sidebar from "./components/sidebar";
import FrameProgressBar from "./components/progressbar";
import {
  ParticleInstances,
  BondInstances,
  SimulationCell,
  PerParticleVectors,
  getCentroid,
} from "./components/particles";
import { Frames, Frame, Player } from "./components/particles";
import { Geometries } from "./components/geometries";
import "./App.css";
import { Plotting } from "./components/plotting";

import { Canvas, useThree } from "@react-three/fiber";
import {
  OrbitControls,
  PerspectiveCamera,
  TrackballControls,
  TransformControls,
} from "@react-three/drei";
import { Button, InputGroup, Form } from "react-bootstrap";
import * as THREE from "three";
import { Line3D, VirtualCanvas } from "./components/lines";
import ControlsBuilder from "./components/transforms";
import { ParticleInfoOverlay, SceneInfoOverlay } from "./components/overlays";
import VectorField from "./components/vectorfield";
import { useColorMode } from "./components/utils";

export default function App() {
  // const [isConnected, setIsConnected] = useState(socket.connected);
  // const [fooEvents, setFooEvents] = useState([]);

  const [selectionSchema, setSelectionSchema] = useState({});
  const [modifierSchema, setModifierSchema] = useState({});
  const [sceneSchema, setSceneSchema] = useState({});
  const [geometrySchema, setGeometrySchema] = useState({});
  const [analysisSchema, setAnalysisSchema] = useState({});
  const [currentFrame, setCurrentFrame] = useState<Frame>({
    arrays: { colors: [], radii: [] },
    calc: null,
    cell: [],
    connectivity: [],
    info: null,
    numbers: [],
    pbc: [],
    positions: [],
    vectors: [],
  });
  const [playing, setPlaying] = useState<boolean>(false);
  const [length, setLength] = useState<number>(0);
  // updated via sockets
  const [step, setStep] = useState<number>(0);
  const [selectedFrames, setSelectedFrames] = useState<Set<number>>(new Set());
  const [selectedIds, setSelectedIds] = useState<Set<number>>(new Set());
  const [bookmarks, setBookmarks] = useState<any>({}); // {name: [step, ...]
  const [points, setPoints] = useState<THREE.Vector3[]>([]);

  const stepFromSocket = useRef<boolean>(false);
  const bookmarksFromSocket = useRef<boolean>(true);
  const selectionFromSocket = useRef<boolean>(true);
  const pointsFromSocket = useRef<boolean>(true);
  const cameraFromSocket = useRef<boolean>(true);

  const [isDrawing, setIsDrawing] = useState<boolean>(false);
  const [selectedPoint, setSelectedPoint] = useState<THREE.Vector3 | null>(
    null,
  );
  const [needsUpdate, setNeedsUpdate] = useState<boolean>(false);
  const [roomName, setRoomName] = useState<string>("");
  const [geometries, setGeometries] = useState<any>([]);
  const [orbitControlsTarget, setOrbitControlsTarget] = useState<THREE.Vector3>(
    new THREE.Vector3(0, 0, 0),
  );
  const [cameraPosition, setCameraPosition] = useState<THREE.Vector3>(
    new THREE.Vector3(0, 0, 0),
  );
  // TODO: initial values are wrong for orbitcontrolstarget and camperaPosition
  // todo give to particles and bonds
  const [colorMode, handleColorMode] = useColorMode();
  const [hoveredId, setHoveredId] = useState<number>(null);
  const [roomConfig, setRoomConfig] = useState({
    arrows: {},
    scene: {},
  });

  const [isAuthenticated, setIsAuthenticated] = useState<boolean>(true);
  const [roomLock, setRoomLock] = useState<boolean>(false);

  // QUEUES
  const [modifierQueue, setModifierQueue] = useState<number>(-1);
  const [selectionQueue, setSelectionQueue] = useState<number>(-1);
  const [geometryQueue, setGeometryQueue] = useState<number>(-1);
  const [analysisQueue, setAnalysisQueue] = useState<number>(-1);

  const [triggerSelection, setTriggerSelection] = useState<boolean>(false);

  const cameraLightRef = useRef<THREE.PointLight>(null);
  const controlsRef = useRef<TransformControls>(null);
  const cameraRef = useRef<THREE.Camera>(null);

  // extension UI elements
  const [tutorialURL, setTutorialURL] = useState<string>("");
  const [showSiMGen, setShowSiMGen] = useState<boolean>(false);

  const [cursorPosition, setCursorPosition] = useState({ x: 0, y: 0 });
  const [lineLength, setLineLength] = useState<number>(0);
  const [showParticleInfo, setShowParticleInfo] = useState<boolean>(false);

  // external useEffects, should be disabled when
  // the input is received via sockets
  sendStep(step, stepFromSocket);
  sendSelection(selectedIds, selectionFromSocket);
  sendBookmarks(bookmarks, bookmarksFromSocket);
  sendCamera(
    {
      position: cameraPosition.toArray(),
      target: orbitControlsTarget.toArray(),
    },
    cameraFromSocket,
  );
  sendPoints(points, pointsFromSocket);

  // if step changes
  useEffect(() => {
    socket.emit("room:frames:get", [step], (frames: Frames) => {
      for (const key in frames) {
        if (frames.hasOwnProperty(key)) {
          const frame: Frame = frames[key]["value"];
          frame.positions = frame.positions.map(
            (position) =>
              new THREE.Vector3(position[0], position[1], position[2]),
          ) as THREE.Vector3[];
        }
        setCurrentFrame(frames[step]["value"]);
        setNeedsUpdate(false); // rename this to something more descriptive
      }
    });
  }, [step, needsUpdate]);

  useEffect(() => {
    // TODO can't be here, because is dependent on the length
    function onFramesRefresh(updatedFrames: number[]) {
      socket.emit("room:length:get", (data: number | string) => {
        // ensure that data is a number
        if (updatedFrames.includes(step)) {
          setNeedsUpdate(true);
        } else if (step >= data) {
          setStep(parseInt(data) - 1);
          // reset selected ids
          setSelectedIds(new Set());
        } else if (roomConfig.scene.frame_update) {
          setStep(parseInt(data) - 1);
        }
        setLength(data);
      });
    }
    socket.on("room:frames:refresh", onFramesRefresh);

    return () => {
      socket.off("room:frames:refresh", onFramesRefresh);
    };
  }, [step, roomConfig.scene]);

  useEffect(() => {
    function onConnect() {
      socket.emit(
        "webclient:connect",
        (data: { name: string; room: string; authenticated: boolean }) => {
          setRoomName(data["room"]);
          setIsAuthenticated(data["authenticated"]);
        },
      );
      console.log("connected");
      // get length
      socket.emit("room:length:get", (data: number) => {
        setLength(data);
        console.log("number of available frames", data);
      });
      // get bookmarks
      socket.emit("room:bookmarks:get", (data: any) => {
        bookmarksFromSocket.current = true;
        setBookmarks(data);
      });
      // // get points
      socket.emit("room:points:get", (data: { [key: string]: number[][] }) => {
        pointsFromSocket.current = true;
        setPoints(data["0"].map((x) => new THREE.Vector3(...x)));
      });
      // get geometries
      socket.emit("room:geometry:get", (data: any) => {
        setGeometries(data);
      });
      // get step
      socket.emit("room:step:get", (data: string) => {
        // this happens only once, we can afford sending the step back.
        // bugfix for missing out modifying the step first in the UI
        // stepFromSocket.current = true;
        setStep(parseInt(data));
      });
      // get selection
      socket.emit("room:selection:get", (data: any) => {
        selectionFromSocket.current = true;
        setSelectedIds(new Set(data[0]));
      });

      // get lock state
      socket.emit("room:lock:get", (data: boolean) => {
        setRoomLock(data);
      });

      // get the config
      socket.emit("room:config:get", (data: any) => {
        setRoomConfig(data);
      });
    }

    function onSelectionSchema(receivedSchema: any) {
      setSelectionSchema(receivedSchema);
    }
    function onModifierSchema(receivedSchema: any) {
      setModifierSchema(receivedSchema);
    }
    function onSceneSchema(receivedSchema: any) {
      setSceneSchema(receivedSchema);
    }
    function onGeometryScheme(receivedSchema: any) {
      setGeometrySchema(receivedSchema);
    }
    function onAnalysisSchema(receivedSchema: any) {
      setAnalysisSchema(receivedSchema);
    }
    // data is {collection_id: [id1, id2, ...]}
    function onRoomSelectionSet(data: any) {
      selectionFromSocket.current = true;
      setSelectedIds(new Set(data[0]));
    }

    function onSetStep(newStep: number) {
      stepFromSocket.current = true;
      setStep(newStep);
    }

    function onGeometries(data: any) {
      setGeometries(data);
    }

    function onBookmarks(data: any) {
      bookmarksFromSocket.current = true;
      setBookmarks(data);
    }

    function onPointsSet(points: { 0: number[][] }) {
      pointsFromSocket.current = true;
      setPoints(points[0].map((point) => new THREE.Vector3(...point)));
    }

    function onModifierQueue(data: number) {
      setModifierQueue(data);
    }

    function onAnalysisQueue(data: number) {
      setAnalysisQueue(data);
    }

    function onGeometryQueue(data: number) {
      setGeometryQueue(data);
    }

    function onSelectionQueue(data: number) {
      setSelectionQueue(data);
    }

    function onModifierRefresh() {
      socket.emit("modifier:schema");
    }
    function onAnalysisRefresh() {
      socket.emit("analysis:schema");
      // scene provides visualisation of vectors, needs to know what is available
      socket.emit("scene:schema");
    }

    function onTutorialURL(data: string) {
      setTutorialURL(data);
    }
    function onShowSiMGen(data: boolean) {
      setShowSiMGen(data);
    }

    function onRoomConfig(data: any) {
      setRoomConfig((prevConfig: any) => ({
        ...prevConfig,
        ...data,
      }));
    }

    function onCameraSet(data: { position: number[]; target: number[] }) {
      cameraFromSocket.current = true;
      setOrbitControlsTarget(new THREE.Vector3(...data.target));
      setCameraPosition(new THREE.Vector3(...data.position));
      if (controlsRef.current && cameraRef.current) {
        controlsRef.current.enabled = false;
        cameraRef.current.position.set(...data.position);
        controlsRef.current.update();
        controlsRef.current.enabled = true;
      }
    }

    function onRoomLockSet(locked: boolean) {
      setRoomLock(locked);
    }

    socket.on("connect", onConnect);
    socket.on("selection:schema", onSelectionSchema);
    socket.on("modifier:schema", onModifierSchema);
    socket.on("scene:schema", onSceneSchema);
    socket.on("geometry:schema", onGeometryScheme);
    socket.on("analysis:schema", onAnalysisSchema);
    socket.on("room:selection:set", onRoomSelectionSet);
    socket.on("room:step:set", onSetStep);
    socket.on("room:geometry:set", onGeometries);
    socket.on("room:bookmarks:set", onBookmarks);
    socket.on("room:modifier:queue", onModifierQueue);
    socket.on("room:analysis:queue", onAnalysisQueue);
    socket.on("room:geometry:queue", onGeometryQueue);
    socket.on("room:selection:queue", onSelectionQueue);
    socket.on("modifier:schema:refresh", onModifierRefresh);
    socket.on("analysis:schema:refresh", onAnalysisRefresh);
    socket.on("room:points:set", onPointsSet);
    socket.on("tutorial:url", onTutorialURL);
    socket.on("showSiMGen", onShowSiMGen);
    socket.on("room:camera:set", onCameraSet);
    socket.on("room:lock:set", onRoomLockSet);
    socket.on("room:config:set", onRoomConfig);

    return () => {
      socket.off("connect", onConnect);
      socket.off("selection:schema", onSelectionSchema);
      socket.off("modifier:schema", onModifierSchema);
      socket.off("scene:schema", onSceneSchema);
      socket.off("geometry:schema", onGeometryScheme);
      socket.off("analysis:schema", onAnalysisSchema);
      socket.off("room:selection:set", onRoomSelectionSet);
      socket.off("room:step:set", onSetStep);
      socket.off("room:geometry:set", onGeometries);
      socket.off("room:bookmarks:set", onBookmarks);
      socket.off("room:modifier:queue", onModifierQueue);
      socket.off("room:analysis:queue", onAnalysisQueue);
      socket.off("room:geometry:queue", onGeometryQueue);
      socket.off("room:selection:queue", onSelectionQueue);
      socket.off("modifier:schema:refresh", onModifierRefresh);
      socket.off("analysis:schema:refresh", onAnalysisRefresh);
      socket.off("room:points:set", onPointsSet);
      socket.off("tutorial:url", onTutorialURL);
      socket.off("showSiMGen", onShowSiMGen);
      socket.off("room:camera:set", onCameraSet);
      socket.off("room:lock:set", onRoomLockSet);
      socket.off("room:config:set", onRoomConfig);
    };
  }, []);

  useEffect(() => {
    // page initialization
    const updateLength = () => {
      socket.emit("room:length:get", (data: number) => {
        setLength(data);
      });
    };

    const handleKeyDown = (event: KeyboardEvent) => {
      // if canvas is not focused, don't do anything
      if (document.activeElement !== document.body) {
        return;
      }
      if (event.key === "ArrowRight") {
        setPlaying(false);
        if (event.shiftKey) {
          // Jump to the bookmark after the current step
          const bookmarkKeys = Object.keys(bookmarks)
            .map(Number)
            .sort((a, b) => a - b);
          const nextBookmark = bookmarkKeys.find((key) => key > step);

          if (nextBookmark !== undefined) {
            setStep(nextBookmark);
          }
        } else {
          if (selectedFrames.size > 0) {
            const nextFrame = Array.from(selectedFrames).find(
              (frame) => frame > step,
            );
            if (nextFrame) {
              setStep(nextFrame);
            } else {
              setStep(Math.min(...selectedFrames));
            }
          } else {
            setStep((prevStep) => (prevStep + 1 < length ? prevStep + 1 : 0));
          }
        }
      } else if (event.key === "ArrowLeft") {
        setPlaying(false);
        if (event.shiftKey) {
          // Jump to the bookmark before the current step
          const bookmarkKeys = Object.keys(bookmarks)
            .map(Number)
            .sort((a, b) => b - a);
          const previousBookmark = bookmarkKeys.find((key) => key < step);

          if (previousBookmark !== undefined) {
            setStep(previousBookmark);
          }
        } else {
          // Move to the previous step, or wrap around to the end
          // check if selectedFrames length is greater than 0, then only jump
          // between selectedFrames
          if (selectedFrames.size > 0) {
            const previousFrame = Array.from(selectedFrames)
              .reverse()
              .find((frame) => frame < step);
            if (previousFrame) {
              setStep(previousFrame);
            } else {
              setStep(Math.max(...selectedFrames));
            }
          } else {
            setStep((prevStep) => (prevStep - 1 >= 0 ? prevStep - 1 : length));
          }
        }
      } else if (event.key == "ArrowUp") {
        // jump 10 percent, or to the end
        setPlaying(false);
        const newStep = Math.min(step + Math.floor(length / 10), length - 1);
        setStep(newStep);
      } else if (event.key == "ArrowDown") {
        // jump 10 percent, or to the beginning
        setPlaying(false);
        const newStep = Math.max(step - Math.floor(length / 10), 0);
        setStep(newStep);
      } else if (event.key == " ") {
        // backspace
        updateLength();
        setPlaying((prev) => !prev);
        if (step == length - 1) {
          setStep(0);
        }
      } else if (event.key == "x") {
        setIsDrawing((prev) => !prev);
      } else if (event.key == "i") {
        setShowParticleInfo((prev) => !prev);
      } else if (event.key == "b") {
        setBookmarks((prev) => {
          const newBookmarks = { ...prev };
          newBookmarks[step] = `Frame ${step}`;
          return newBookmarks;
        });
      } else if (event.key == "a") {
        if (event.ctrlKey) {
          setSelectedIds(
            new Set([
              ...Array.from(
                { length: currentFrame.positions.length },
                (_, i) => i,
              ),
            ]),
          );
        }
      } else if (event.key == "Backspace" || event.key == "Delete") {
        // check if shift is pressed
        if (event.shiftKey) {
          if (selectedPoint !== null) {
            const newPoints = points.filter(
              (point) => point.distanceTo(selectedPoint) > 0.1,
            );
            setSelectedPoint(null);
            setPoints(newPoints);
          } else if (points.length > 0) {
            // pop last point from points
            setSelectedPoint(null);
            setPoints(points.slice(0, points.length - 1));
          }
        } else {
          if (selectedIds.size > 0) {
            socket.emit("modifier:run", {
              method: { discriminator: "Delete" },
            });
          }
        }
      } else if (event.key == "c") {
        if (selectedPoint !== null) {
          setOrbitControlsTarget(selectedPoint);
        } else {
          if (currentFrame.positions.length > 0) {
            const center = getCentroid(currentFrame.positions, selectedIds);
            setOrbitControlsTarget(center);
          }
        }
      } else if (event.key == "o") {
        const origin = {
          position: [-10, -10, -10],
          target: getCentroid(currentFrame.positions, new Set()),
        };
        setOrbitControlsTarget(new THREE.Vector3(...origin.target));
        setCameraPosition(new THREE.Vector3(...origin.position));
        if (controlsRef.current && cameraRef.current) {
          controlsRef.current.enabled = false;
          cameraRef.current.position.set(...origin.position);
          controlsRef.current.update();
          controlsRef.current.enabled = true;
        }
      }
    };

    // Add the event listener
    window.addEventListener("keydown", handleKeyDown);

    // Clean up the event listener on unmount
    return () => {
      window.removeEventListener("keydown", handleKeyDown);
    };
  }, [
    length,
    step,
    points,
    selectedPoint,
    bookmarks,
    currentFrame,
    selectedIds,
  ]);

  useEffect(() => {
    const updateCursorPosition = (event) => {
      setCursorPosition({ x: event.clientX, y: event.clientY });
    };

    window.addEventListener("mousemove", updateCursorPosition);

    return () => {
      window.removeEventListener("mousemove", updateCursorPosition);
    };
  }, []);

  const onDragOver = (event) => {
    event.preventDefault();
  };

  const onDrop = async (event) => {
    event.preventDefault();

    const file = event.dataTransfer.files[0];
    if (!file) {
      console.error("No file was dropped");
      return;
    }

    try {
      const arrayBuffer = await file.arrayBuffer();
      const content = new Uint8Array(arrayBuffer);

      // send the file to the server
      socket.emit("room:upload:file", {
        content: Array.from(content),
        filename: file.name,
      });
    } catch (error) {
      console.error("Error reading file:", error);
    }
  };

  const onPointerMissed = () => {
    setSelectedPoint(null);
    setSelectedIds(new Set());
  };

  const setSceneSettings = (data: any) => {
    setRoomConfig((prev) => ({
      ...prev,
      scene: data,
    }));
    socket.emit("room:config:set", { scene: data });
  };

  return (
    <>
      <div className="canvas-container" onDragOver={onDragOver} onDrop={onDrop}>
        {roomConfig.scene.controls !== undefined && (
          <Canvas onPointerMissed={onPointerMissed}>
            <PerspectiveCamera
              ref={cameraRef}
              makeDefault
              near={roomConfig["scene"]["camera_near"]}
              far={roomConfig["scene"]["camera_far"]}
            />
            <pointLight
              ref={cameraLightRef}
              position={[0, 0, 0]}
              decay={0}
              intensity={Math.PI}
              castShadow
            />
            {roomConfig["scene"]["vectorfield"] &&
              currentFrame.vectors !== undefined && (
                <VectorField
                  vectors={currentFrame.vectors}
                  arrowsConfig={{
                    rescale: roomConfig["scene"].vector_scale,
                    ...roomConfig.arrows,
                  }}
                />
              )}
            <ParticleInstances
              frame={currentFrame}
              selectedIds={selectedIds}
              setSelectedIds={setSelectedIds}
              isDrawing={isDrawing}
              points={points}
              setPoints={setPoints}
              setOrbitControlsTarget={setOrbitControlsTarget}
              hoveredId={hoveredId}
              setHoveredId={setHoveredId}
              setTriggerSelection={setTriggerSelection}
              sceneSettings={roomConfig["scene"]}
            />
            <BondInstances
              frame={currentFrame}
              selectedIds={selectedIds}
              hoveredId={hoveredId}
              sceneSettings={roomConfig["scene"]}
            />
            {roomConfig["scene"]["simulation_box"] && (
              <SimulationCell frame={currentFrame} colorMode={colorMode} />
            )}
            {roomConfig["scene"].controls === "OrbitControls" && (
              <OrbitControls
                ref={controlsRef}
                enableDamping={false}
                target={orbitControlsTarget}
                onChange={(e) => {
                  if (!e) return;
                  const camera = e.target.object;
                  if (cameraLightRef.current) {
                    cameraLightRef.current.position
                      .copy(camera.position)
                      .sub(orbitControlsTarget)
                      .normalize()
                      .add(camera.position);
                  }
                  setCameraPosition(camera.position);
                }}
                makeDefault
              />
            )}
            {roomConfig["scene"].controls === "TrackballControls" && (
              <TrackballControls
                ref={controlsRef}
                target={orbitControlsTarget}
                staticMoving={true}
                onChange={(e) => {
                  if (!e) return;
                  const camera = e.target.object;
                  if (cameraLightRef.current) {
                    cameraLightRef.current.position
                      .copy(camera.position)
                      .sub(orbitControlsTarget)
                      .normalize()
                      .add(camera.position);
                  }
                  setCameraPosition(camera.position);
                }}
                makeDefault
              />
            )}
            <Player
              playing={playing}
              togglePlaying={setPlaying}
              step={step}
              setStep={setStep}
              fps={roomConfig["scene"].fps}
              loop={roomConfig["scene"]["Animation Loop"]}
              length={length}
              selectedFrames={selectedFrames}
            />
            <Line3D
              points={points}
              setPoints={setPoints}
              setSelectedPoint={setSelectedPoint}
              isDrawing={isDrawing}
              colorMode={colorMode}
              hoveredId={hoveredId}
              setIsDrawing={setIsDrawing}
              setLineLength={setLineLength}
            />
            <ControlsBuilder
              points={points}
              setPoints={setPoints}
              selectedPoint={selectedPoint}
              setSelectedPoint={setSelectedPoint}
            />
            <Geometries
              geometries={geometries}
              isDrawing={isDrawing}
              setHoveredId={setHoveredId}
              setPoints={setPoints}
            />
            <VirtualCanvas
              setPoints={setPoints}
              isDrawing={isDrawing}
              points={points}
              hoveredId={hoveredId}
              setHoveredId={setHoveredId}
            />
            {roomConfig["scene"].vectors != "" && (
              <PerParticleVectors
                frame={currentFrame}
                property={roomConfig["scene"].vectors}
                colorMode={colorMode}
                arrowsConfig={{
                  rescale: roomConfig["scene"].vector_scale,
                  ...roomConfig.arrows,
                }}
              ></PerParticleVectors>
            )}
          </Canvas>
        )}
      </div>
      <div className="App">
        <HeadBar
          room={roomName}
          colorMode={colorMode}
          handleColorMode={handleColorMode}
          setIsDrawing={setIsDrawing}
          setGeometries={setGeometries}
          setPoints={setPoints}
          isDrawing={isDrawing}
          tutorialURL={tutorialURL}
          showSiMGen={showSiMGen}
          modifierQueue={modifierQueue}
          isAuthenticated={isAuthenticated}
          roomLock={roomLock}
        />
        <Sidebar
          selectionSchema={selectionSchema}
          modifierSchema={modifierSchema}
          sceneSchema={sceneSchema}
          geometrySchema={geometrySchema}
          analysisSchema={analysisSchema}
          sceneSettings={roomConfig["scene"]}
          setSceneSettings={setSceneSettings}
          modifierQueue={modifierQueue}
          selectionQueue={selectionQueue}
          geometryQueue={geometryQueue}
          analysisQueue={analysisQueue}
          triggerSelection={triggerSelection}
          setTriggerSelection={setTriggerSelection}
          colorMode={colorMode}
          setStep={setStep}
        />
        <FrameProgressBar
          length={length}
          step={step}
          setStep={setStep}
          bookmarks={bookmarks}
          setBookmarks={setBookmarks}
          selectedFrames={selectedFrames}
          setSelectedFrames={setSelectedFrames}
        />
        <Plotting setStep={setStep} setSelectedFrames={setSelectedFrames} />
        {showParticleInfo && (
          <>
            <ParticleInfoOverlay
              show={hoveredId !== null || isDrawing}
              info={{
                ...(hoveredId !== null && {
                  "Particle ID": hoveredId,
                  "Atomic Number": currentFrame.numbers[hoveredId],
                }),
                ...(isDrawing && { Line: `${lineLength.toFixed(2)} Å` }),
              }}
              position={cursorPosition}
            />
            <SceneInfoOverlay
              frame={currentFrame}
              setShowParticleInfo={setShowParticleInfo}
            />
          </>
        )}
      </div>
    </>
  );
}
