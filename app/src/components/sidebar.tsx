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
  storageKey: string;
}

const SidebarMenu: React.FC<SidebarMenuProps> = ({
  schema,
  onSubmit,
  queuePosition,
  trigger,
  setTrigger,
  visible,
  storageKey,
}) => {
  const editorRef = useRef<HTMLDivElement>(null);
  const editor = useRef<JSONEditor>(null);
  const [userInput, setUserInput] = useState<any>(null);

  function submitEditor() {
    if (onSubmit) {
      onSubmit(editor.current.getValue());
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

  useEffect(() => {
    if (Object.keys(schema).length === 0) {
      return;
    }
    let setNewEditor = false;
    if (!editor.current) {
      setNewEditor = true;
    }
    if (editor.current) {
      // console.log("editor exists and ready: ", editor.current.ready);
      setNewEditor = !editor.current.ready;
    }

    if (setNewEditor && editorRef.current) {
      editor.current = new JSONEditor(editorRef.current, {
        schema: schema,
      });
      setNewEditor = false;
      let triggerChange = false;

      editor.current.on("change", () => {
        if (editor.current.ready && triggerChange) {
          const editorValue = editor.current.getValue();
          editor.current.validate();
          setUserInput(editorValue);
          // localStorage.setItem(storageKey, JSON.stringify(editorValue));
        }
      });

      editor.current.on("ready", () => {
        // const userInput = localStorage.getItem(storageKey);
        // if (userInput) {
        if (userInput) {
          editor.current.setValue(userInput);
        }
        // if (Object.keys(JSON.parse(userInput)).length > 0) {
        //   editor.current.setValue(JSON.parse(userInput));
        // }
        // }
        setTimeout(() => {
          triggerChange = true;
        }, 1000);
      });
    }
    return () => {
      if (editor.current) {
        editor.current.destroy();
      }
    };
  }, [schema]);

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

function AnalysisMenu({
  selectionSchema,
  onSubmit,
  plotData,
  setPlotData,
  queuePosition,
  visible,
  storageKey,
}: {
  selectionSchema: any;
  onSubmit: any;
  plotData: any;
  setPlotData: any;
  queuePosition: number;
  visible: boolean;
  storageKey: string;
}) {
  const editorRef = useRef<HTMLDivElement>(null);
  const figureRef = useRef<HTMLDivElement>(null);
  // the values that are currently in the editor, only used because
  // we don't have a reference to the editor object
  const [currentEditorValue, setCurrentEditorValue] = useState<any>({});
  // TODO: these are lost when the component is unmounted

  useEffect(() => {
    const handleFigure = (data) => {
      try {
        const parsedData = JSON.parse(data);
        setPlotData(parsedData);
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
      onSubmit(currentEditorValue);
    }
  }

  useEffect(() => {
    if (editorRef.current) {
      const editor = new JSONEditor(editorRef.current, {
        schema: selectionSchema,
        theme: "bootstrap5",
        iconlib: "fontawesome5",
        show_errors: "always",
      });
      editor.on("change", () => {
        editor.validate();
        setCurrentEditorValue(editor.getValue());
        // if (onChange) {
        //   onChange(editor.getValue());
        // }
      });

      return () => {
        if (editorRef.current) {
          editorRef.current.innerHTML = "";
        }
      };
    }
  }, [selectionSchema]); // this does only trigger if selectionSchema changes

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
          <Button onClick={submitEditor} disabled={queuePosition >= 0}>
            {queuePosition > 0 && `Queue position: ${queuePosition}`}
            {queuePosition == 0 && `Running`}
            {queuePosition < 0 && `Submit`}
          </Button>
        </Card.Body>
        <Card.Footer>
          {plotData.data.length > 0 && (
            <Plot
              ref={figureRef}
              data={plotData.data}
              layout={plotData.layout}
              frames={plotData.frames}
              config={plotData.config}
              style={{ width: "100%", height: "100%" }}
            />
          )}
        </Card.Footer>
      </Card>
    </>
  );
}

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
}) {
  const [visibleOption, setVisibleOption] = useState<string>("");
  const [plotData, setPlotData] = useState({
    data: [],
    layout: {},
    frames: [],
    config: {},
  });

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
        storageKey="selection-menu"
      />
      <SidebarMenu
        schema={modifierSchema}
        onSubmit={(data: any) => {
          socket.emit("modifier:run", data);
        }}
        queuePosition={modifierQueue}
        visible={visibleOption == "interaction"}
        storageKey="interaction-menu"
      />
      <SidebarMenu
        schema={sceneSchema}
        onSubmit={setSceneSettings}
        queuePosition={-1}
        visible={visibleOption == "scene"}
        storageKey="scene-menu"
      />
      <SidebarMenu
        schema={geometrySchema}
        onSubmit={(data: any) => {
          socket.emit("geometry:run", data);
        }}
        queuePosition={geometryQueue}
        visible={visibleOption == "geometry"}
        storageKey="geometry-menu"
      />
      <AnalysisMenu
        selectionSchema={analysisSchema}
        onSubmit={(data: any) => {
          socket.emit("analysis:run", data);
        }}
        plotData={plotData}
        setPlotData={setPlotData}
        queuePosition={analysisQueue}
        visible={visibleOption == "analysis"}
        storageKey="analysis-menu"
      />
    </>
  );
}

export default SideBar;
