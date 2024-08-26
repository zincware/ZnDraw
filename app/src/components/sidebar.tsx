import { Button, Navbar, Nav, Card } from "react-bootstrap";
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

JSONEditor.defaults.options.theme = "bootstrap5";
JSONEditor.defaults.options.iconlib = "fontawesome5";
JSONEditor.defaults.options.object_background = "bg-body-color";
JSONEditor.defaults.options.disable_edit_json = true;
JSONEditor.defaults.options.disable_properties = true;
JSONEditor.defaults.options.disable_collapse = true;
JSONEditor.defaults.options.no_additional_properties = true;
JSONEditor.defaults.options.keep_oneof_values = false;
JSONEditor.defaults.editors.object.options.titleHidden = true;

interface SidebarMenuProps {
  schema: any;
  onSubmit: any;
  queuePosition: number;
  trigger?: boolean; // Mark trigger as optional
  setTrigger?: (value: boolean) => void; // Mark setTrigger as optional
  visible: boolean;
  useSubmit?: boolean; // provide a submit button or trigger on change
  closeMenu?: () => void;
}

const useJSONEditor = (
  schema: any,
  setUserInput: (value: any) => void,
  useSubmit: boolean,
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
      let created_trigger = false;

      // on ready, validate and set user input
      JSONEditorRef.current.on("ready", () => {
        if (useSubmit) {
          // when using the submit button, we need to set the user input on ready
          // otherwise, it could be None.
          console.log("Setting user input on ready");
          if (JSONEditorRef.current.validate()) {
            const editorValue = JSONEditorRef.current.getValue();
            setUserInput(editorValue);
          }
        }
      });

      JSONEditorRef.current.on("change", () => {
        if (JSONEditorRef.current.ready) {
          if (created_trigger) {
            if (JSONEditorRef.current.validate()) {
              const editorValue = JSONEditorRef.current.getValue();
              setUserInput(editorValue);
            }
          } else {
            // skip first trigger
            created_trigger = true;
          }
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
  useSubmit,
  closeMenu,
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

  const editorRef = useJSONEditor(schema, setUserInput, useSubmit);

  useEffect(() => {
    if (!useSubmit && userInput !== null) {
      submitEditor();
    }
  }, [userInput, useSubmit]);

  return (
    <Card
      className="rounded-0 border-start-0 overflow-y-auto rounded-end"
      style={{
        position: "fixed",
        top: "50px",
        left: "50px",
        height: "100%",
        maxWidth: "35%",
        margin: 0,
        padding: 0,
        display: visible ? "block" : "none",
      }}
    >
      <Card.Header
        className="d-flex justify-content-between align-items-center"
        style={{
          backgroundColor: "inherit", // Use the same background color as the rest of the card
        }}
      >
        <Card.Title>{schema.title}</Card.Title>
        <Button variant="close" className="ms-auto" onClick={closeMenu} />
      </Card.Header>
      <Card.Body style={{ marginTop: -30, paddingBottom: 80 }}>
        <div ref={editorRef}></div>
        {useSubmit && (
          <Button onClick={submitEditor} disabled={queuePosition >= 0}>
            {queuePosition > 0 && `Queue position: ${queuePosition}`}
            {queuePosition == 0 && `Running`}
            {queuePosition < 0 && `Submit`}
          </Button>
        )}
      </Card.Body>
    </Card>
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
  setStep,
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
  setStep: any;
}) {
  const [visibleOption, setVisibleOption] = useState<string>("");
  useEffect(() => {
    if (visibleOption !== "") {
      // emit visibleOption:schema e.g. selection:schema
      socket.emit(`${visibleOption}:schema`);
    }
  }, [visibleOption]);

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
                    visibleOption != "modifier" ? "modifier" : "",
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
        useSubmit={true}
        closeMenu={() => setVisibleOption("")}
      />
      <SidebarMenu
        schema={modifierSchema}
        onSubmit={(data: any) => {
          socket.emit("modifier:run", data);
        }}
        queuePosition={modifierQueue}
        visible={visibleOption == "modifier"}
        useSubmit={true}
        closeMenu={() => setVisibleOption("")}
      />
      <SidebarMenu
        schema={sceneSchema}
        onSubmit={setSceneSettings}
        queuePosition={-1}
        visible={visibleOption == "scene"}
        useSubmit={false}
        closeMenu={() => setVisibleOption("")}
      />
      <SidebarMenu
        schema={geometrySchema}
        onSubmit={(data: any) => {
          socket.emit("geometry:run", data);
        }}
        queuePosition={geometryQueue}
        visible={visibleOption == "geometry"}
        useSubmit={true}
        closeMenu={() => setVisibleOption("")}
      />
      <SidebarMenu
        schema={analysisSchema}
        onSubmit={(data: any) => {
          socket.emit("analysis:run", data);
        }}
        queuePosition={analysisQueue}
        visible={visibleOption == "analysis"}
        useSubmit={true}
        closeMenu={() => setVisibleOption("")}
      />
    </>
  );
}

export default SideBar;
