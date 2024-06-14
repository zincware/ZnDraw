import { Button, Navbar, Nav, Card, Collapse } from "react-bootstrap";
import { BtnTooltip } from "./tooltips";
import {
  FaRegChartBar,
  FaCube,
  FaRegHandPointer,
  FaRegMap,
  FaGithub,
} from "react-icons/fa";
import { FaCircleNodes } from "react-icons/fa6";
import { useState, useRef } from "react";
import { socket } from "../socket";
import { useEffect } from "react";

import { JSONEditor } from "@json-editor/json-editor";
import Plot from "react-plotly.js";
import { Rnd } from "react-rnd";

JSONEditor.defaults.options.theme = "bootstrap5";
JSONEditor.defaults.options.iconlib = "fontawesome5";
JSONEditor.defaults.options.object_background = "bg-body-color";
JSONEditor.defaults.options.disable_edit_json = true;
JSONEditor.defaults.options.disable_properties = true;
JSONEditor.defaults.options.disable_collapse = true;
JSONEditor.defaults.options.no_additional_properties = true;
JSONEditor.defaults.options.keep_oneof_values = false;

interface SidebarMenuProps {
  schema: any;
  onSubmit: any;
  queuePosition: number;
  trigger?: boolean; // Mark trigger as optional
  setTrigger?: (value: boolean) => void; // Mark setTrigger as optional
  visible: boolean;
}

interface AnalysisMenuProps extends SidebarMenuProps {
  showPlotsCard: boolean;
  setPlotData: (data: any) => void;
  colorMode: string;
  setShowPlotsCard: (value: boolean) => void;
}

const useJSONEditor = (
  schema: any,
  userInput: any,
  setUserInput: (value: any) => void,
) => {
  const editorRef = useRef<HTMLDivElement>(null);
  const JSONEditorRef = useRef<JSONEditor | null>(null);

  useEffect(() => {
    if (Object.keys(schema).length === 0) {
      return;
    }

    if (editorRef.current) {
      JSONEditorRef.current = new JSONEditor(editorRef.current, {
        schema: schema,
      });

      JSONEditorRef.current.on("change", () => {
        if (JSONEditorRef.current.ready) {
          const editorValue = JSONEditorRef.current.getValue();
          JSONEditorRef.current.validate();
          setUserInput(editorValue);
        }
      });

      JSONEditorRef.current.on("ready", () => {
        if (userInput) {
          JSONEditorRef.current.setValue(userInput);
        }
      });
    }
    return () => {
      if (JSONEditorRef.current) {
        JSONEditorRef.current.destroy();
      }
    };
  }, [schema]);

  return editorRef;
};

const SidebarMenu: React.FC<SidebarMenuProps> = ({
  schema,
  onSubmit,
  queuePosition,
  trigger,
  setTrigger,
  visible,
}) => {
  const [userInput, setUserInput] = useState<any>(null);

  function submitEditor() {
    if (onSubmit) {
      onSubmit(userInput);
    }
  }

  useEffect(() => {
    if (trigger) {
      submitEditor();
      if (setTrigger) {
        setTrigger(false);
      }
    }
  }, [trigger]);

  const editorRef = useJSONEditor(schema, userInput, setUserInput);

  return (
    <Card
      className="rounded-0 border-start-0 overflow-y-auto rounded-end"
      style={{
        position: "fixed",
        top: "50px",
        left: "50px",
        height: "100%",
        display: visible ? "block" : "none",
      }}
    >
      <Card.Body className="pt-0 pb-5 my-0">
        <div ref={editorRef}></div>
        <Button onClick={submitEditor} disabled={queuePosition >= 0}>
          {queuePosition > 0 && `Queue position: ${queuePosition}`}
          {queuePosition == 0 && `Running`}
          {queuePosition < 0 && `Submit`}
        </Button>
      </Card.Body>
    </Card>
  );
};

const AnalysisMenu: React.FC<AnalysisMenuProps> = ({
  schema,
  onSubmit,
  setPlotData,
  queuePosition,
  visible,
  colorMode,
  showPlotsCard,
  setShowPlotsCard,
}) => {
  const [userInput, setUserInput] = useState<any>(null);
  const [plotStyle, setPlotStyle] = useState<any>({
    width: "100%",
    height: "100%",
  });
  const editorRef = useJSONEditor(schema, userInput, setUserInput);

  useEffect(() => {
    const handleFigure = (data) => {
      try {
        const parsedData = JSON.parse(data);
        setPlotData(parsedData);
        setShowPlotsCard(true);
      } catch (error) {
        console.error("Error parsing JSON data: ", error);
      }
    };

    socket.on("analysis:figure:set", handleFigure);

    return () => {
      socket.off("analysis:figure:set", handleFigure);
    };
  }, []);

  function submitEditor() {
    if (onSubmit) {
      onSubmit(userInput);
    }
  }

  return (
    <>
      <Card
        className="rounded-0 border-start-0 overflow-y-auto rounded-end"
        style={{
          position: "fixed",
          top: "50px",
          left: "50px",
          height: "100%",
          width: "50%",
          display: visible ? "block" : "none",
        }}
      >
        <Card.Body className="pt-0 pb-5 my-0">
          <div ref={editorRef}></div>
          <Button onClick={submitEditor} disabled={queuePosition >= 0} className="mx-2">
            {queuePosition > 0 && `Queue position: ${queuePosition}`}
            {queuePosition == 0 && `Running`}
            {queuePosition < 0 && `Submit`}
          </Button>
          <Button
            onClick={() => {
              setShowPlotsCard(!showPlotsCard);
            }}>
            {showPlotsCard ? "Hide Figure" : "Show Figure"}
            </Button>
        </Card.Body>
      </Card>
    </>
  );
};


const PlotsCard = ({
  plotData,
  setPlotData,
  colorMode,
  showPlotsCard,
  setShowPlotsCard,
}: {
  plotData: any;
  setPlotData: any;
  colorMode: string;
  showPlotsCard: boolean;
  setShowPlotsCard: any;
}) => {
  const [plotStyle, setPlotStyle] = useState<any>({
    width: "100%",
    height: "100%",
  });
  const [renderKey, setRenderKey] = useState<number>(0);
  const cardRef = useRef<any>(null);

  useEffect(() => {
    if (plotData) {
      const newPlotData = { ...plotData };
      newPlotData.layout.paper_bgcolor = colorMode === "dark" ? "rgba(0,0,0, 0)" : "rgba(255,255,255, 0)";
      setPlotData(newPlotData);
    }
    const newPlotStyle = { ...plotStyle };
    newPlotStyle.filter = colorMode === "dark" ? "invert(75%) hue-rotate(180deg)" : "";
    setPlotStyle(newPlotStyle);
  }, [colorMode]);

  const onResize = () => {
    if (cardRef.current) {
      const newPlotData = { ...plotData };
      newPlotData.layout.width = cardRef.current.clientWidth;
      newPlotData.layout.height = cardRef.current.clientHeight - 50;
      setPlotData(newPlotData);
      setRenderKey((prevKey) => prevKey + 1);
    }
  };

  return (
    <Rnd
      default={{
        x: 100,
        y: -100,
        width: 400,
        height: 400,
      }}
      style={{ zIndex: 1000, padding: 0, margin: 0, display: showPlotsCard ? "block" : "none" }}
      onResize={onResize}
    >
      <Card
        style={{
          margin: 0,
          padding: 0,
          width: "100%",
          height: "100%",
        }}
        ref={cardRef}
      >
        <Card.Header className="d-flex justify-content-between align-items-center">
          <Card.Title>Analysis Figure</Card.Title>
          <Button variant="close" onClick={() => setShowPlotsCard(false)} />
        </Card.Header>
        <Card.Body>
          {plotData.data.length > 0 && (
            <Plot
              key={renderKey}
              data={plotData.data}
              layout={plotData.layout}
              frames={plotData.frames}
              config={plotData.config}
              style={plotStyle}
            />
          )}
        </Card.Body>
      </Card>
    </Rnd>
  );
};

function SideBar({
  selectionSchema,
  modifierSchema,
  sceneSchema,
  geometrySchema,
  analysisSchema,
  sceneSettings,
  setSceneSettings,
  modifierQueue,
  selectionQueue,
  analysisQueue,
  geometryQueue,
  triggerSelection,
  setTriggerSelection,
  colorMode,
}: {
  selectionSchema: any;
  modifierSchema: any;
  sceneSchema: any;
  geometrySchema: any;
  analysisSchema: any;
  sceneSettings: any;
  setSceneSettings: any;
  modifierQueue: number;
  selectionQueue: number;
  analysisQueue: number;
  geometryQueue: number;
  triggerSelection: boolean;
  setTriggerSelection: any;
  colorMode: string;
}) {
  const [visibleOption, setVisibleOption] = useState<string>("");
  const [plotData, setPlotData] = useState({
    data: [],
    layout: {},
    frames: [],
    config: {},
  });
  const [showPlotsCard, setShowPlotsCard] = useState<boolean>(false);

  // if any menu is open and you click escape, close it
  useEffect(() => {
    function handleKeyDown(event: KeyboardEvent) {
      if (event.key === "Escape") {
        setVisibleOption("");
      }
    }

    document.addEventListener("keydown", handleKeyDown);

    return () => {
      document.removeEventListener("keydown", handleKeyDown);
    };
  }, []);

  return (
    <>
      <Navbar className="flex-column bg-body-primary bg-tertiary">
        <Card
          border="tertiary py-0 px-0"
          className="rounded-0 border-top-0"
          style={{
            position: "fixed",
            top: "50px",
            left: "0",
            height: "100%",
            width: "50px",
          }}
        >
          <BtnTooltip text="Selection" placement="right">
            <Nav className="mx-auto my-1">
              <Button
                variant="outline-tertiary"
                onClick={() =>
                  setVisibleOption(
                    visibleOption != "selection" ? "selection" : "",
                  )
                }
              >
                <FaRegHandPointer />
              </Button>
            </Nav>
          </BtnTooltip>
          <BtnTooltip text="Interaction" placement="right">
            <Nav className="mx-auto my-1">
              <Button
                variant="outline-tertiary"
                onClick={() =>
                  setVisibleOption(
                    visibleOption != "interaction" ? "interaction" : "",
                  )
                }
              >
                <FaCircleNodes />
              </Button>
            </Nav>
          </BtnTooltip>
          <BtnTooltip text="Scene" placement="right">
            <Nav className="mx-auto my-1">
              <Button
                variant="outline-tertiary"
                onClick={() =>
                  setVisibleOption(visibleOption != "scene" ? "scene" : "")
                }
              >
                <FaCube />
              </Button>
            </Nav>
          </BtnTooltip>
          <BtnTooltip text="Geometry" placement="right">
            <Nav className="mx-auto my-1">
              <Button
                variant="outline-tertiary"
                onClick={() =>
                  setVisibleOption(
                    visibleOption != "geometry" ? "geometry" : "",
                  )
                }
              >
                <FaRegMap />
              </Button>
            </Nav>
          </BtnTooltip>
          <BtnTooltip text="Analysis" placement="right">
            <Nav className="mx-auto my-1">
              <Button
                variant="outline-tertiary"
                onClick={() =>
                  setVisibleOption(
                    visibleOption != "analysis" ? "analysis" : "",
                  )
                }
              >
                <FaRegChartBar />
              </Button>
            </Nav>
          </BtnTooltip>
          <BtnTooltip text="View on GitHub" placement="right">
            <Nav className="mx-auto my-1">
              <Button
                variant="outline-tertiary"
                href="https://github.com/zincware/ZnDraw"
                target="_blank"
              >
                <FaGithub />
              </Button>
            </Nav>
          </BtnTooltip>
        </Card>
      </Navbar>
      <SidebarMenu
        schema={selectionSchema}
        onSubmit={(data: any) => {
          socket.emit("selection:run", data);
        }}
        queuePosition={selectionQueue}
        trigger={triggerSelection}
        setTrigger={setTriggerSelection}
        visible={visibleOption == "selection"}
      />
      <SidebarMenu
        schema={modifierSchema}
        onSubmit={(data: any) => {
          socket.emit("modifier:run", data);
        }}
        queuePosition={modifierQueue}
        visible={visibleOption == "interaction"}
      />
      <SidebarMenu
        schema={sceneSchema}
        onSubmit={setSceneSettings}
        queuePosition={-1}
        visible={visibleOption == "scene"}
      />
      <SidebarMenu
        schema={geometrySchema}
        onSubmit={(data: any) => {
          socket.emit("geometry:run", data);
        }}
        queuePosition={geometryQueue}
        visible={visibleOption == "geometry"}
      />
      <AnalysisMenu
        schema={analysisSchema}
        onSubmit={(data: any) => {
          socket.emit("analysis:run", data);
        }}
        setPlotData={setPlotData}
        queuePosition={analysisQueue}
        visible={visibleOption == "analysis"}
        colorMode={colorMode}
        showPlotsCard={showPlotsCard}
        setShowPlotsCard={setShowPlotsCard}
      />
      <PlotsCard
        setPlotData={setPlotData}
        plotData={plotData}
        colorMode={colorMode}
        showPlotsCard={showPlotsCard}
        setShowPlotsCard={setShowPlotsCard}
      />
    </>
  );
}

export default SideBar;
