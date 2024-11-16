import React, { useState, useEffect, useRef } from "react";
import { socket, client } from "./socket";
import { Pathtracer } from "@react-three/gpu-pathtracer";
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
  setupConfig,
} from "./components/api";
import HeadBar from "./components/headbar";
import Sidebar from "./components/sidebar";
import FrameProgressBar from "./components/progressbar";
import {
  ParticleInstances,
  BondInstances,
  SimulationCell,
  PerParticleVectors,
} from "./components/particles";
import { getCentroid } from "./components/particlesEditor";
import { Frames, Frame, Player } from "./components/particles";
import { Geometries } from "./components/geometries";
import "./App.css";
import { Plotting } from "./components/plotting";
import * as znsocket from "znsocket";

import { Canvas, useFrame } from "@react-three/fiber";
import CameraAndControls from "./components/cameraAndControls";
import {
  OrbitControls,
  PerspectiveCamera,
  OrthographicCamera,
  TrackballControls,
  TransformControls,
  Box,
  CameraControls,
  Environment,
  Sphere,
} from "@react-three/drei";
import { Button, InputGroup, Form } from "react-bootstrap";
import * as THREE from "three";
import { Line3D, VirtualCanvas } from "./components/lines";
import ControlsBuilder from "./components/transforms";
import { ParticleInfoOverlay, SceneInfoOverlay } from "./components/overlays";
import VectorField from "./components/vectorfield";
import { useColorMode } from "./components/utils";
import { IndicesState } from "./components/utils";
import { Floor } from "./components/floor";

export default function App() {
  // const [isConnected, setIsConnected] = useState(socket.connected);
  // const [fooEvents, setFooEvents] = useState([]);

  const [connected, setConnected] = useState<boolean>(false);
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
  const [cameraAndControls, setCameraAndControls] = useState<any>({
    camera: new THREE.Vector3(0, 0, 0),
    target: new THREE.Vector3(0, 0, 0),
  });
  // TODO: initial values are wrong for orbitcontrolstarget and camperaPosition
  // todo give to particles and bonds
  const [colorMode, handleColorMode] = useColorMode();
  const [hoveredId, setHoveredId] = useState<number>(-1);
  // UPDATE THESE using `vis.config` on the Python side
  const [roomConfig, setRoomConfig] = useState({
    arrows: {
      colormap: [
        [-0.5, 0.9, 0.5],
        [0.0, 0.9, 0.5],
      ],
      normalize: true,
      colorrange: [0, 1.0],
      scale_vector_thickness: false,
      opacity: 1.0,
    },
    scene: {
      fps: 30,
      material: "MeshStandardMaterial",
      particle_size: 1.0,
      bond_size: 1.0,
      animation_loop: false,
      simulation_box: true,
      vectorfield: true,
      controls: "OrbitControls",
      vectors: [],
      vector_scale: 1.0,
      selection_color: "#ffa500",
      camera: "PerspectiveCamera",
      camera_near: 0.1,
      camera_far: 300,
      frame_update: true,
      crosshair: false,
      floor: false,
      synchronize_camera: true,
    },
    PathTracer: {
      enabled: false,
      environment: "studio",
      metalness: 0.7,
      roughness: 0.2,
      clearcoat: 0.0,
      clearcoatRoughness: 0.0,
    },
  });

  const [isAuthenticated, setIsAuthenticated] = useState<boolean>(true);
  const [roomLock, setRoomLock] = useState<boolean>(false);

  // QUEUES
  // TODO: fix
  const [modifierQueue, setModifierQueue] = useState<number>(-1);

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
  setupConfig(token, setRoomConfig);
  setupBookmarks(token, setBookmarks, bookmarks);
  setupPoints(token, setPoints, points);
  setupSelection(token, setSelectedIds, selectedIds);
  setupStep(token, setStep, step);
  setupCamera(
    token,
    cameraAndControls,
    setCameraAndControls,
    roomConfig["scene"]["synchronize_camera"],
  );
  setupFrames(
    token,
    step,
    setCurrentFrame,
    currentFrame,
    setLength,
    setStep,
    roomConfig["scene"]["frame_update"],
  );
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

      socket.emit("room:token:get", (data: string) => {
        setToken(data);
      });
    }

    function onDisconnect() {
      setConnected(false);
    }

    function onTutorialURL(data: string) {
      setTutorialURL(data);
    }
    function onShowSiMGen(data: boolean) {
      setShowSiMGen(data);
    }

    function onRoomLockSet(locked: boolean) {
      setRoomLock(locked);
    }

    socket.on("connect", onConnect);
    socket.on("disconnect", onDisconnect);
    socket.on("tutorial:url", onTutorialURL);
    socket.on("showSiMGen", onShowSiMGen);
    socket.on("room:lock:set", onRoomLockSet);

    return () => {
      socket.off("connect", onConnect);
      socket.off("disconnect", onDisconnect);
      socket.off("tutorial:url", onTutorialURL);
      socket.off("showSiMGen", onShowSiMGen);
      socket.off("room:lock:set", onRoomLockSet);
    };
  }, []);

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
            const queue = new znsocket.Dict({
              client: client,
              key: "queue:" + token + ":" + "modifier",
            });
            queue["Delete"] = {};
            socket.emit("room:worker:run");
          }
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

  // reduce selection, if selected points is reduced
  useEffect(() => {
    if (selectedIds.size > 0) {
      const newSelectedIds = new Set(
        Array.from(selectedIds).filter(
          (id) => id < currentFrame.positions.length,
        ),
      );
      // if the selection is reduced, update the selection
      if (newSelectedIds.size < selectedIds.size) {
        setSelectedIds(newSelectedIds);
      }
    }
  }, [currentFrame, selectedIds]);

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

  return (
    <>
      <div className="canvas-container" onDragOver={onDragOver} onDrop={onDrop}>
        <Canvas onPointerMissed={onPointerMissed} shadows>
          <CameraAndControls
            roomConfig={roomConfig}
            cameraAndControls={cameraAndControls}
            setCameraAndControls={setCameraAndControls}
            currentFrame={currentFrame}
            selectedIds={selectedIds}
            colorMode={colorMode}
          />
          <Pathtracer enabled={roomConfig.PathTracer.enabled}>
            {roomConfig.PathTracer.enabled &&
              roomConfig.PathTracer.environment !== "none" && (
                <Environment preset={roomConfig.PathTracer.environment} />
              )}

            {roomConfig["scene"].floor ? (
              <>
                <Floor colorMode={colorMode} roomConfig={roomConfig} />
                <directionalLight
                  position={[0, 100, 0]}
                  intensity={1.0}
                  castShadow
                  shadow-mapSize-width={roomConfig["scene"]["camera_far"] * 10} // Adjust the width of the shadow map
                  shadow-mapSize-height={roomConfig["scene"]["camera_far"] * 10} // Adjust the height of the shadow map
                  shadow-camera-near={10} // Adjust the near clipping plane of the shadow camera
                  shadow-camera-far={800} // Adjust the far clipping plane of the shadow camera
                  shadow-camera-left={-1 * roomConfig["scene"]["camera_far"]} // Set the left boundary for the shadow camera frustum
                  shadow-camera-right={roomConfig["scene"]["camera_far"]} // Set the right boundary for the shadow camera frustum
                  shadow-camera-top={roomConfig["scene"]["camera_far"]} // Set the top boundary for the shadow camera frustum
                  shadow-camera-bottom={-1 * roomConfig["scene"]["camera_far"]} // Set the bottom boundary for the shadow camera frustum
                />
              </>
            ) : (
              <directionalLight position={[0, 100, 0]} intensity={1.0} />
            )}

            {roomConfig["scene"]["vectorfield"] &&
              currentFrame.vectors !== undefined && (
                <VectorField
                  vectors={currentFrame.vectors}
                  pathTracingSettings={roomConfig.PathTracer}
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
              setPoints={setPoints}
              setHoveredId={setHoveredId}
              sceneSettings={roomConfig["scene"]}
              token={token}
              highlight=""
              visibleIndices={undefined}
              setFrame={setCurrentFrame}
              pathTracingSettings={roomConfig.PathTracer}
            />
            {!roomConfig.PathTracer.enabled && (
              <>
                <ParticleInstances
                  frame={currentFrame}
                  selectedIds={selectedIds}
                  setSelectedIds={setSelectedIds}
                  isDrawing={isDrawing}
                  setPoints={setPoints}
                  setHoveredId={setHoveredId}
                  sceneSettings={roomConfig["scene"]}
                  token={token}
                  visibleIndices={hoveredId}
                  highlight={"backside"}
                  setFrame={setCurrentFrame}
                />
                <ParticleInstances
                  frame={currentFrame}
                  selectedIds={selectedIds}
                  setSelectedIds={setSelectedIds}
                  isDrawing={isDrawing}
                  setPoints={setPoints}
                  setHoveredId={setHoveredId}
                  sceneSettings={roomConfig["scene"]}
                  token={token}
                  visibleIndices={selectedIds}
                  highlight={"selection"}
                  setFrame={setCurrentFrame}
                />
                <ParticleInstances
                  frame={currentFrame}
                  selectedIds={selectedIds}
                  setSelectedIds={setSelectedIds}
                  isDrawing={isDrawing}
                  setPoints={setPoints}
                  setHoveredId={setHoveredId}
                  sceneSettings={roomConfig["scene"]}
                  token={token}
                  visibleIndices={
                    new Set(currentFrame.constraints?.[0]?.indices)
                  }
                  highlight={"constraint"}
                  setFrame={setCurrentFrame}
                />
                <BondInstances
                  frame={currentFrame}
                  visibleIndices={selectedIds}
                  highlight="selection"
                  sceneSettings={roomConfig["scene"]}
                />
              </>
            )}
            <BondInstances
              frame={currentFrame}
              visibleIndices={undefined}
              highlight=""
              sceneSettings={roomConfig["scene"]}
              pathTracingSettings={roomConfig.PathTracer}
            />
            {roomConfig["scene"]["simulation_box"] &&
              !roomConfig.PathTracer.enabled && (
                <SimulationCell frame={currentFrame} colorMode={colorMode} />
              )}
            <Player
              playing={playing}
              togglePlaying={setPlaying}
              step={step}
              setStep={setStep}
              fps={roomConfig["scene"].fps}
              loop={roomConfig["scene"]["animation_loop"]}
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
            {roomConfig["scene"].vectors[0] &&
              roomConfig["scene"].vectors.map((vector) => (
                <PerParticleVectors
                  frame={currentFrame}
                  property={vector}
                  colorMode={colorMode}
                  arrowsConfig={{
                    rescale: roomConfig["scene"].vector_scale,
                    ...roomConfig.arrows,
                  }}
                  pathTracingSettings={roomConfig.PathTracer}
                  key={vector}
                ></PerParticleVectors>
              ))}
          </Pathtracer>
        </Canvas>
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
        <Sidebar token={token} />
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
