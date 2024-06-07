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
} from "./components/particles";
import { Frames, Frame, Player } from "./components/particles";
import { Geometries } from "./components/geometries";
import "./App.css";

import { Canvas, useThree } from "@react-three/fiber";
import {
  OrbitControls,
  PerspectiveCamera,
  TransformControls,
} from "@react-three/drei";
import { Button, InputGroup, Form } from "react-bootstrap";
import * as THREE from "three";
import { Line3D, VirtualCanvas } from "./components/lines";
import ControlsBuilder from "./components/transforms";

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
  });
  const [playing, setPlaying] = useState<boolean>(false);
  const [length, setLength] = useState<number>(0);
  const [sceneSettings, setSceneSettings] = useState({
    fps: 60,
    "Animation Loop": false,
    simulation_box: false,
  });
  // updated via sockets
  const [step, setStep] = useState<number>(0);
  const [selectedIds, setSelectedIds] = useState<Set<number>>(new Set());
  const [bookmarks, setBookmarks] = useState<any>({}); // {name: [step, ...]
  const [points, setPoints] = useState<THREE.Vector3[]>([]);

  const stepFromSocket = useRef<boolean>(true); // true to avoid first render trigger
  const bookmarksFromSocket = useRef<boolean>(true);
  const selectionFromSocket = useRef<boolean>(true);
  const pointsFromSocket = useRef<boolean>(true); 

  const [isDrawing, setIsDrawing] = useState<boolean>(false);
  const [selectedPoint, setSelectedPoint] = useState<THREE.Vector3 | null>(
    null,
  );
  const [needsUpdate, setNeedsUpdate] = useState<boolean>(false);
  const [roomName, setRoomName] = useState<string>("");
  const [geometries, setGeometries] = useState<any>([]);
  const [orbitControlsTarget, setOrbitControlsTarget] = useState<
    [number, number, number]
  >([0, 0, 0]);
  const [cameraPosition, setCameraPosition] = useState<
    [number, number, number]
  >([0, 0, 0]);
  // TODO: initial values are wrong for orbitcontrolstarget and camperaPosition
  // todo give to particles and bonds
  const [colorMode, setColorMode] = useState<string>("light");
  const [hoveredId, setHoveredId] = useState(null);

  // QUEUES
  const [modifierQueue, setModifierQueue] = useState<number>(-1);
  const [selectionQueue, setSelectionQueue] = useState<number>(-1);
  const [geometryQueue, setGeometryQueue] = useState<number>(-1);
  const [analysisQueue, setAnalysisQueue] = useState<number>(-1);

  const cameraLightRef = useRef<THREE.PointLight>(null);
  const cameraRef = useRef<THREE.Camera>(null);

  // external useEffects, should be disabled when
  // the input is received via sockets
  sendStep(step, stepFromSocket);
  sendSelection(selectedIds, selectionFromSocket);
  sendBookmarks(bookmarks, bookmarksFromSocket);
  sendCamera({ position: cameraPosition, target: orbitControlsTarget });
  sendPoints(points, pointsFromSocket);

  // if step changes
  useEffect(() => {
    socket.emit("room:frames:get", [step], (data: Frames) => {
      const key = Object.keys(data)[0];
      if (key) {
        setCurrentFrame(data[key]);
        setNeedsUpdate(false);
        // TODO: cache
      }
    });
  }, [step, needsUpdate]);

  useEffect(() => {
    // TODO can't be here, because is dependent on the length
    function onFramesRefresh(updatedFrames: number[]) {
      socket.emit("room:length:get", (data: number) => {
        if (updatedFrames.includes(step)) {
          console.log("step is in updated frames", step);
          setNeedsUpdate(true);
        } else if (step >= data) {
          console.log("step is out of bounds", step);
          setStep(data - 1);
          // reset selected ids
          setSelectedIds(new Set());
        }
        setLength(data);
      });
    }
    socket.on("room:frames:refresh", onFramesRefresh);

    return () => {
      socket.off("room:frames:refresh", onFramesRefresh);
    };
  }, [step]);

  useEffect(() => {
    function onConnect() {
      socket.emit(
        "webclient:connect",
        (data: { name: string; room: string }) => {
          setRoomName(data["room"]);
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
      socket.emit("room:points:get", (data: any) => {
        pointsFromSocket.current = true;
        setPoints(data["0"]);
      });
      // get geometries
      socket.emit("room:geometry:get", (data: any) => {
        setGeometries(data);
      });
      // get step
      socket.emit("room:step:get", (data: number) => {
        stepFromSocket.current = true;
        setStep(data);
      }); 
      // get selection
      socket.emit("room:selection:get", (data: any) => {
        selectionFromSocket.current = true;
        setSelectedIds(new Set(data));
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
      console.log("geometries", data);
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
    };
  }, []);

  useEffect(() => {
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
      if (event.key == "ArrowRight") {
        if (event.shiftKey) {
          // jump to the bookmark after the current step
          const keys = Object.keys(bookmarks).map((x) => parseInt(x));
          keys.sort((a, b) => a - b);
          const next_bookmark = keys.find((x) => x > step);
          if (next_bookmark) {
            setStep(next_bookmark);
          }
        } else if (step + 1 < length) {
          setStep((prev) => prev + 1);
        } else {
          setStep(0);
        }
      } else if (event.key == "ArrowLeft") {
        if (event.shiftKey) {
          // jump to the bookmark before the current step
          const keys = Object.keys(bookmarks).map((x) => parseInt(x));
          keys.sort((a, b) => b - a);
          const next_bookmark = keys.find((x) => x < step);
          if (next_bookmark) {
            setStep(next_bookmark);
          }
        } else {
          if (step - 1 >= 0) {
            setStep((prev) => prev - 1);
          } else {
            setStep(length - 1);
          }
        }
      } else if (event.key == "ArrowUp") {
        // jump 10 percent, or to the end
        const newStep = Math.min(step + Math.floor(length / 10), length - 1);
        setStep(newStep);
      } else if (event.key == "ArrowDown") {
        // jump 10 percent, or to the beginning
        const newStep = Math.max(step - Math.floor(length / 10), 0);
        setStep(newStep);
      } else if (event.key == " ") {
        updateLength();
        setPlaying((prev) => !prev);
        if (step == length - 1) {
          setStep(0);
        }
      } else if (event.key == "x") {
        setIsDrawing((prev) => !prev);
      } else if (event.key == "b") {
        setBookmarks((prev) => {
          const newBookmarks = { ...prev };
          newBookmarks[step] = `Frame ${step}`;
          return newBookmarks;
        });
      }
    };

    // Add the event listener
    window.addEventListener("keydown", handleKeyDown);

    // Clean up the event listener on unmount
    return () => {
      window.removeEventListener("keydown", handleKeyDown);
    };
  }, [length, step]);

  return (
    <>
      <div className="canvas-container">
        <Canvas onPointerMissed={() => setSelectedPoint(null)}>
          <PerspectiveCamera ref={cameraRef} />
          {/* <ambientLight intensity={Math.PI / 20}/> */}
          <pointLight
            ref={cameraLightRef}
            position={[2, 2, -10]} // camera position is [0, 0, 5]
            decay={0}
            intensity={Math.PI}
            castShadow
          />

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
          />
          <BondInstances
            frame={currentFrame}
            selectedIds={selectedIds}
            setSelectedIds={setSelectedIds}
            isDrawing={isDrawing}
            points={points}
            setPoints={setPoints}
          />
          {sceneSettings["simulation_box"] && (
            <SimulationCell frame={currentFrame} colorMode={colorMode} />
          )}
          <OrbitControls
            enableDamping={false}
            target={orbitControlsTarget}
            onChange={(e) => {
              if (!e) return;
              const camera = e.target.object;
              if (cameraLightRef.current) {
                cameraLightRef.current.position.set(2, 2, 0);
                cameraLightRef.current.position.add(camera.position);
              }
              setCameraPosition(camera.position.toArray());
            }}
            makeDefault
          />
          <Player
            playing={playing}
            togglePlaying={setPlaying}
            step={step}
            setStep={setStep}
            fps={sceneSettings.fps}
            loop={sceneSettings["Animation Loop"]}
            length={length}
          />
          <Line3D
            points={points}
            setPoints={setPoints}
            setSelectedPoint={setSelectedPoint}
            isDrawing={isDrawing}
            colorMode={colorMode}
            hoveredId={hoveredId}
            setIsDrawing={setIsDrawing}
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
          />
        </Canvas>
      </div>
      <div className="App">
        <HeadBar
          room={roomName}
          colorMode={colorMode}
          setColorMode={setColorMode}
          setIsDrawing={setIsDrawing}
          setGeometries={setGeometries}
          setPoints={setPoints}
        />
        <Sidebar
          selectionSchema={selectionSchema}
          modifierSchema={modifierSchema}
          sceneSchema={sceneSchema}
          geometrySchema={geometrySchema}
          analysisSchema={analysisSchema}
          sceneSettings={sceneSettings}
          setSceneSettings={setSceneSettings}
          modifierQueue={modifierQueue}
          selectionQueue={selectionQueue}
          geometryQueue={geometryQueue}
          analysisQueue={analysisQueue}
        />
        <FrameProgressBar
          length={length}
          step={step}
          setStep={setStep}
          bookmarks={bookmarks}
          setBookmarks={setBookmarks}
        />
      </div>
    </>
  );
}
