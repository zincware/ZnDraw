import { SetStateAction, useEffect, useState, useMemo } from "react";
import {
  Navbar,
  Nav,
  Container,
  Button,
  Modal,
  Card,
  ToggleButton,
} from "react-bootstrap";
import {
  FaCode,
  FaDownload,
  FaHandSparkles,
  FaMoon,
  FaRegClipboard,
  FaSun,
  FaTerminal,
  FaUpload,
} from "react-icons/fa";
import Markdown from "react-markdown";
import { Rnd } from "react-rnd";
import { BtnTooltip } from "./tooltips";
import { socket } from "../socket";

import { FaArrowRotateRight, FaCircleInfo, FaPencil } from "react-icons/fa6";

function ConsoleWindow({ text }: { text: string[] }) {
  return (
    <Rnd
      default={{
        x: window.innerWidth / 2 - 210,
        y: window.innerHeight / 2 - 200,
        width: 200,
        height: "100px",
      }}
      maxHeight={"150px"}
      minWidth={"200px"}
      style={{ zIndex: 1000, padding: 0, margin: 0 }}
    >
      <Card
        style={{
          margin: 0,
          padding: 0,
          // background: "rgba(255, 255, 255, 0.85)",
          // backdropFilter: "blur(5px)",
        }}
      >
        <Card.Header>
          <Card.Title>Console</Card.Title>
        </Card.Header>
        <Card.Body>
          <div style={{ overflowY: "auto", height: "100px" }}>
            {text.map((line, idx) => (
              <p key={idx}>{line}</p>
            ))}
          </div>
        </Card.Body>
      </Card>
    </Rnd>
  );
}

function HelpModel(props: any) {
  const helpMD = `
- **play / pause**: \`keypress space\`
- **frame forwards / backwards**: \`keypress ▶\\◀\`
- **jump forwards / backwards**: \`keypress ▲\\▼\`
- **center camera around selected particle**: \`keypress C\`
- **select multiple particles**: \`keydown shift\`
- **show particle index**: \`keydown I\`
- **toggle canvas movement**: \`keypress F\`
- **duplicate selected anchor point**: \`keypress D\`
- **switch transform control mode**: \`keypress T\`
- **select / scale drawing geometry**: \`keypress shift + mousewheel\`
- **add bookmark at current step**: \`keypress B\`
- **jump between bookmarks**: \`shift + keypress ▶\\◀\`
- **remove bookmark**: \`shift + mouse click\`
- **toggle drawing mode**: \`keypress X\`
`;

  return (
    <Modal {...props} aria-labelledby="contained-modal-title-vcenter" size="lg">
      <Modal.Header closeButton>
        <Modal.Title id="contained-modal-title-vcenter">
          ZnDraw Help
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>
        <Container>
          <Markdown>{helpMD}</Markdown>
        </Container>
      </Modal.Body>
      <Modal.Footer>
        <Button onClick={props.onHide}>Close</Button>
      </Modal.Footer>
    </Modal>
  );
}

function ConnectModal({ show, onHide, room }) {
  // const url = window.location.href.replace(/\/$/, "");
  // for testing the socketio url is different and hard coded
  const serverUrl = window.location.origin;
  // const serverUrl = "http://localhost:1235";

  const pythonCode = `from zndraw import ZnDraw

vis = ZnDraw(url="${serverUrl}", token="${room}")`;

  const helpMD = `
You can connect a Python kernel to this ZnDraw session using

\`\`\`python
${pythonCode}
\`\`\`
  `;

  const copyToClipboard = () => {
    navigator.clipboard.writeText(pythonCode);
    onHide();
  };

  return (
    <Modal
      show={show}
      onHide={onHide}
      aria-labelledby="contained-modal-title-vcenter"
      size="lg"
    >
      <Modal.Header closeButton>
        <Modal.Title id="contained-modal-title-vcenter">
          ZnDraw Help
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>
        <Container>
          <Markdown>{helpMD}</Markdown>
        </Container>
      </Modal.Body>
      <Modal.Footer>
        <Button onClick={copyToClipboard}>
          <FaRegClipboard /> Copy to Clipboard
        </Button>
        <Button onClick={onHide}>Close</Button>
      </Modal.Footer>
    </Modal>
  );
}

function RefreshModal({ show, onHide, room }) {
  const urlWithRoom =
    window.location.href.replace(/\/$/, "") + `/token/${room}`;

  return (
    <Modal
      show={show}
      onHide={onHide}
      aria-labelledby="contained-modal-title-vcenter"
      size="lg"
    >
      <Modal.Header closeButton>
        <Modal.Title id="contained-modal-title-vcenter">
          Reload Scene
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>
        <Container>
          <p>
            You are about the reload the scene. The current scene will be
            available as long as some client is connected, otherwise the data
            will be deleted. You can access the scene by visiting <br />
            <a href={urlWithRoom}>{urlWithRoom}</a>
          </p>
          <Button
            variant="outline-primary"
            onClick={() => navigator.clipboard.writeText(urlWithRoom)}
          >
            copy URL to clipboard
          </Button>
        </Container>
      </Modal.Body>
      <Modal.Footer>
        <Button onClick={onHide}>Cancel</Button>
        <Button href={`${window.location.href}reset`}>Create new Scene</Button>
      </Modal.Footer>
    </Modal>
  );
}

const HeadBar = ({
  room,
  colorMode,
  setColorMode,
  setIsDrawing,
  setGeometries,
  setPoints,
  isDrawing,
}: {
  room: string;
  colorMode: string;
  setColorMode: any;
  setIsDrawing: any;
  setGeometries: any;
  setPoints: any;
  isDrawing: boolean;
}) => {
  const [helpModalShow, setHelpModalShow] = useState(false);
  const [connectModalShow, setConnectModalShow] = useState(false);
  const [refreshModalShow, setRefreshModalShow] = useState(false);
  const [consoleShow, setConsoleShow] = useState(false);
  const [consoleText, setConsoleText] = useState<string[]>([]);

  useEffect(() => {
    const handleConsoleText = (text: string) => {
      const date = new Date().toISOString();
      const msg = date.slice(0, -8) + " " + text;
      setConsoleText((prev) => [msg, ...prev]);
    };
    socket.on("room:log", handleConsoleText);
    return () => {
      socket.off("room:log", handleConsoleText);
    };
  }, []);

  const handleColorMode = () => {
    setColorMode(colorMode === "light" ? "dark" : "light");
    document.documentElement.setAttribute(
      "data-bs-theme",
      colorMode === "light" ? "dark" : "light",
    );
  };

  const handleRemovePointsGeometries = () => {
    console.log("remove points and geometries");
    socket.emit("room:geometry:set", []);
    socket.emit("room:point:set", []);
    setGeometries([]);
    setPoints([]);
  };

  return (
    <>
      <Navbar
        expand="md"
        className="bg-body-tertiary fixed-top"
        style={{ height: 50 }}
      >
        <Container fluid>
          <Navbar.Brand>ZnDraw</Navbar.Brand>
          <Navbar.Toggle aria-controls="basic-navbar-nav" />
          <Navbar.Collapse id="basic-navbar-nav">
            <Nav className="me-auto">
              <BtnTooltip text="Reset Scene">
                <Button
                  variant="outline-danger"
                  className="mx-1"
                  onClick={() => {
                    setRefreshModalShow((prev: boolean) => !prev);
                  }}
                >
                  <FaArrowRotateRight />
                </Button>
              </BtnTooltip>
              <BtnTooltip text="Activate Drawing Tool">
                <ToggleButton
                  variant="outline-primary"
                  className="mx-1"
                  value="1"
                  id="toggle-drawing"
                  active={isDrawing}
                  onClick={(e) => {
                    setIsDrawing((prev: boolean) => !prev);
                  }}
                >
                  <FaPencil />
                </ToggleButton>
              </BtnTooltip>
              <BtnTooltip text="Remove all guiding points and geometries">
                <Button
                  variant="outline-primary"
                  className="mx-1"
                  onClick={handleRemovePointsGeometries}
                >
                  <FaHandSparkles />
                </Button>
              </BtnTooltip>
            </Nav>
            <Nav className="ms-auto">
              <BtnTooltip text="Access console">
                <ToggleButton
                  variant="outline-primary"
                  className="mx-1"
                  value="1"
                  id="toggle-console"
                  active={consoleShow}
                  onClick={() => {
                    setConsoleShow((prev: boolean) => !prev);
                  }}
                >
                  <FaTerminal />
                </ToggleButton>
              </BtnTooltip>
              <BtnTooltip text="Python access">
                <Button
                  variant="outline-primary"
                  className="mx-1"
                  onClick={() => setConnectModalShow(true)}
                >
                  <FaCode />
                </Button>
              </BtnTooltip>
              {/* <BtnTooltip text="Upload file (max 1 MB)">
                <Button variant="outline-primary" className="mx-1">
                  <FaUpload />
                </Button>
              </BtnTooltip>
              <BtnTooltip text="Download scene">
                <Button variant="outline-primary" className="mx-1">
                  <FaDownload />
                </Button>
              </BtnTooltip> */}
              <BtnTooltip text="Help">
                <Button
                  variant="outline-primary"
                  className="mx-1"
                  onClick={() => setHelpModalShow(true)}
                >
                  <FaCircleInfo />
                </Button>
              </BtnTooltip>
              <BtnTooltip text="Switch Colormode">
                <Button
                  variant="outline-primary"
                  className="mx-1"
                  onClick={handleColorMode}
                >
                  {colorMode === "light" ? <FaSun /> : <FaMoon />}
                </Button>
              </BtnTooltip>
              {/* <Button variant="outline-primary" className="mx-1">
                <FaUsers />
              </Button> */}
            </Nav>
          </Navbar.Collapse>
        </Container>
      </Navbar>
      <HelpModel show={helpModalShow} onHide={() => setHelpModalShow(false)} />
      <ConnectModal
        show={connectModalShow}
        onHide={() => setConnectModalShow(false)}
        room={room}
      />
      <RefreshModal
        show={refreshModalShow}
        onHide={() => setRefreshModalShow(false)}
        room={room}
      />
      {consoleShow && <ConsoleWindow text={consoleText} />}
    </>
  );
};

export default HeadBar;