import React, { useState, useEffect, useRef } from "react";
import { socket, client } from "./socket";
import {
  setupBookmarks,
  setupPoints,
  setupSelection,
  setupStep,
  setupCamera,
  setupFrames,
  setupFigures,
  setupGeometries,
  setupMessages,
} from "./components/api";
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

import { Canvas, useThree, useFrame } from "@react-three/fiber";
import {
  OrbitControls,
  PerspectiveCamera,
  OrthographicCamera,
  TrackballControls,
  TransformControls,
  Box,
  CameraControls,
} from "@react-three/drei";
import { Button, InputGroup, Form } from "react-bootstrap";
import * as THREE from "three";
import { Line3D, VirtualCanvas } from "./components/lines";
import ControlsBuilder from "./components/transforms";
import { ParticleInfoOverlay, SceneInfoOverlay } from "./components/overlays";
import VectorField from "./components/vectorfield";
import { useColorMode } from "./components/utils";
import { IndicesState } from "./components/utils";

const MoveCameraTarget = ({
  controlsRef,
  colorMode,
}: {
  controlsRef: any;
  colorMode: string;
}) => {
  const controlsCrosshairRef = useRef<THREE.Object3D>(null);
  const shortDimension = 0.05;
  const longDimension = 0.5;

  // Update the controlsCrosshair position to match the orbit controls target
  useFrame(() => {
    if (controlsCrosshairRef.current && controlsRef.current) {
      const crosshair = controlsCrosshairRef.current;
      const target = controlsRef.current.target;
      crosshair.position.copy(target);
    }
  });

  return (
    <group ref={controlsCrosshairRef}>
      {/* X axis box */}
      <Box scale={[longDimension, shortDimension, shortDimension]}>
        <meshStandardMaterial
          color={colorMode == "light" ? "#454b66" : "#f5fdc6"}
        />
      </Box>

      {/* Y axis box */}
      <Box scale={[shortDimension, longDimension, shortDimension]}>
        <meshStandardMaterial
          color={colorMode == "light" ? "#454b66" : "#f5fdc6"}
        />
      </Box>

      {/* Z axis box */}
      <Box scale={[shortDimension, shortDimension, longDimension]}>
        <meshStandardMaterial
          color={colorMode == "light" ? "#454b66" : "#f5fdc6"}
        />
      </Box>
    </group>
  );
};

export default function App() {
  // const [isConnected, setIsConnected] = useState(socket.connected);
  // const [fooEvents, setFooEvents] = useState([]);

  const [connected, setConnected] = useState<boolean>(false);
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
  const [selectedFrames, setSelectedFrames] = useState<IndicesState>({
    active: true,
    indices: new Set<number>(),
  });
  const [selectedIds, setSelectedIds] = useState<Set<number>>(new Set());
  const [bookmarks, setBookmarks] = useState<any>({}); // {name: [step, ...]
  const [points, setPoints] = useState<THREE.Vector3[]>([]);

  const [isDrawing, setIsDrawing] = useState<boolean>(false);
  const [selectedPoint, setSelectedPoint] = useState<THREE.Vector3 | null>(
    null,
  );
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

  const [cameraRoll, setCameraRoll] = useState<number>(0); // or undefined

  // extension UI elements
  const [tutorialURL, setTutorialURL] = useState<string>("");
  const [showSiMGen, setShowSiMGen] = useState<boolean>(false);

  const [cursorPosition, setCursorPosition] = useState({ x: 0, y: 0 });
  const [lineLength, setLineLength] = useState<number>(0);
  const [showParticleInfo, setShowParticleInfo] = useState<boolean>(false);
  const [addPlotsWindow, setAddPlotsWindow] = useState<number>(0); // make this bool!
  const [updatedPlotsList, setUpdatedPlotsList] = useState<string[]>([]);
  const [messages, setMessages] = useState<string[]>([]);

  const [token, setToken] = useState<string>("");
  setupBookmarks(token, setBookmarks, bookmarks);
  setupPoints(token, setPoints, points);
  setupSelection(token, setSelectedIds, selectedIds);
  setupStep(token, setStep, step);
  setupCamera(
    token,
    cameraPosition,
    orbitControlsTarget,
    setCameraPosition,
    setOrbitControlsTarget,
    controlsRef,
    cameraRef,
  );
  setupFrames(token, step, setCurrentFrame, setLength, setStep);
  setupFigures(token, setUpdatedPlotsList);
  setupGeometries(token, setGeometries, geometries);
  setupMessages(token, setMessages, messages);

  useEffect(() => {
    function onConnect() {
      setConnected(true);
      socket.emit(
        "webclient:connect",
        (data: { name: string; room: string; authenticated: boolean }) => {
          setRoomName(data["room"]);
          setIsAuthenticated(data["authenticated"]);
        },
      );
      console.log("connected");

      // get lock state
      socket.emit("room:lock:get", (data: boolean) => {
        setRoomLock(data);
      });

      // get the config
      socket.emit("room:config:get", (data: any) => {
        setRoomConfig(data);
      });

      socket.emit("room:token:get", (data: string) => {
        setToken(data);
      });
    }

    function onDisconnect() {
      setConnected(false);
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

    function onRoomLockSet(locked: boolean) {
      setRoomLock(locked);
    }

    socket.on("connect", onConnect);
    socket.on("disconnect", onDisconnect);
    socket.on("selection:schema", onSelectionSchema);
    socket.on("modifier:schema", onModifierSchema);
    socket.on("scene:schema", onSceneSchema);
    socket.on("geometry:schema", onGeometryScheme);
    socket.on("analysis:schema", onAnalysisSchema);
    socket.on("room:modifier:queue", onModifierQueue);
    socket.on("room:analysis:queue", onAnalysisQueue);
    socket.on("room:geometry:queue", onGeometryQueue);
    socket.on("room:selection:queue", onSelectionQueue);
    socket.on("modifier:schema:refresh", onModifierRefresh);
    socket.on("analysis:schema:refresh", onAnalysisRefresh);
    socket.on("tutorial:url", onTutorialURL);
    socket.on("showSiMGen", onShowSiMGen);
    socket.on("room:lock:set", onRoomLockSet);
    socket.on("room:config:set", onRoomConfig);

    return () => {
      socket.off("connect", onConnect);
      socket.off("disconnect", onDisconnect);
      socket.off("selection:schema", onSelectionSchema);
      socket.off("modifier:schema", onModifierSchema);
      socket.off("scene:schema", onSceneSchema);
      socket.off("geometry:schema", onGeometryScheme);
      socket.off("analysis:schema", onAnalysisSchema);
      socket.off("room:modifier:queue", onModifierQueue);
      socket.off("room:analysis:queue", onAnalysisQueue);
      socket.off("room:geometry:queue", onGeometryQueue);
      socket.off("room:selection:queue", onSelectionQueue);
      socket.off("modifier:schema:refresh", onModifierRefresh);
      socket.off("analysis:schema:refresh", onAnalysisRefresh);
      socket.off("tutorial:url", onTutorialURL);
      socket.off("showSiMGen", onShowSiMGen);
      socket.off("room:lock:set", onRoomLockSet);
      socket.off("room:config:set", onRoomConfig);
    };
  }, []);

  // camera roll

  useEffect(() => {
    if (controlsRef.current && cameraRef.current) {
      const camera = cameraRef.current;
      if (camera) {
        controlsRef.current.enabled = false;
        // y direction
        var yDir = new THREE.Vector3(0, 1, 0);
        if (cameraRoll === null) {
          camera.up.copy(yDir);
        } else {
          // test case to roll the camera normal to screen

          // direction camera is looking to
          var looksTo = new THREE.Vector3();
          camera.getWorldDirection(looksTo);

          // direction perpendicular to both yDir and looksTo
          var b = new THREE.Vector3();
          b.crossVectors(yDir, looksTo).normalize();

          // direction perpendicular to both looksTo and b
          var n = new THREE.Vector3();
          n.crossVectors(looksTo, b).normalize();

          // make a circle in the plane with vectors b and n
          n.multiplyScalar(Math.cos(cameraRoll)).add(
            b.multiplyScalar(Math.sin(cameraRoll)),
          );

          // set camera up
          camera.up.set(n.x, n.y, n.z);
        }

        controlsRef.current.update();
        controlsRef.current.enabled = true;
        cameraRef.current.updateProjectionMatrix();
      }
    }
  }, [cameraRoll]);

  useEffect(() => {
    // page initialization

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
          if (selectedFrames.indices.size > 0 && selectedFrames.active) {
            const nextFrame = Array.from(selectedFrames.indices).find(
              (frame) => frame > step,
            );
            if (nextFrame) {
              setStep(nextFrame);
            } else {
              setStep(Math.min(...selectedFrames.indices));
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
          if (selectedFrames.indices.size > 0 && selectedFrames.active) {
            const previousFrame = Array.from(selectedFrames.indices)
              .reverse()
              .find((frame) => frame < step);
            if (previousFrame) {
              setStep(previousFrame);
            } else {
              setStep(Math.max(...selectedFrames.indices));
            }
          } else {
            setStep((prevStep) =>
              prevStep - 1 >= 0 ? prevStep - 1 : length - 1,
            );
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
        // updateLength();
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
          // cameraRef.current.up.set(0, 1, 0); // original up vector
          setCameraRoll(null);
          controlsRef.current.update();
          controlsRef.current.enabled = true;
          cameraRef.current.updateProjectionMatrix();
        }
      } else if (event.key == "r") {
        const roll = Math.PI / 100;
        if (event.ctrlKey) {
          setCameraRoll((prev) => prev - roll);
        } else {
          setCameraRoll((prev) => prev + roll);
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
            {roomConfig["scene"].camera === "PerspectiveCamera" && (
              <PerspectiveCamera
                ref={cameraRef}
                makeDefault
                near={roomConfig["scene"]["camera_near"]}
                far={roomConfig["scene"]["camera_far"]}
                position={[-10, -10, -10]}
              />
            )}
            {roomConfig["scene"].camera === "OrthographicCamera" && (
              <OrthographicCamera
                ref={cameraRef}
                makeDefault
                near={roomConfig["scene"]["camera_near"]}
                far={roomConfig["scene"]["camera_far"]}
                position={[-10, -10, -10]}
                zoom={10}
              />
            )}
            <pointLight
              ref={cameraLightRef}
              position={[-11, -11, -11]}
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
                  setCameraPosition(new THREE.Vector3().copy(camera.position));
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
                  setCameraPosition(new THREE.Vector3().copy(camera.position));
                }}
                makeDefault
              />
            )}
            {roomConfig["scene"].crosshair && (
              <MoveCameraTarget
                controlsRef={controlsRef}
                colorMode={colorMode}
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
          setAddPlotsWindow={setAddPlotsWindow}
          messages={messages}
          setMessages={setMessages}
          token={token}
          step={step}
          selection={selectedIds}
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
          connected={connected}
        />
        <Plotting
          setStep={setStep}
          setSelectedFrames={setSelectedFrames}
          addPlotsWindow={addPlotsWindow}
          setSelectedIds={setSelectedIds}
          step={step}
          updatedPlotsList={updatedPlotsList}
          token={token}
        />
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
