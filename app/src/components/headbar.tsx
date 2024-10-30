import React, {
  SetStateAction,
  useEffect,
  useState,
  useMemo,
  useRef,
  forwardRef,
} from "react";
import Select from 'react-select';
import {
  Navbar,
  Nav,
  Container,
  Button,
  Modal,
  Card,
  ToggleButton,
  InputGroup,
  Form,
} from "react-bootstrap";
import remarkGfm from "remark-gfm";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import "katex/dist/katex.min.css"; // `rehype-katex` does not import the CSS for you
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import {
  oneDark,
  oneLight,
} from "react-syntax-highlighter/dist/esm/styles/prism";
import {
  FaCode,
  FaDownload,
  FaFilm,
  FaHandSparkles,
  FaLock,
  FaLockOpen,
  FaMoon,
  FaRegClipboard,
  FaRocket,
  FaSun,
  FaTerminal,
  FaUpload,
  FaPlus,
} from "react-icons/fa";
import Markdown from "react-markdown";
import { Rnd } from "react-rnd";
import { BtnTooltip } from "./tooltips";
import { socket, client } from "../socket";
import * as znsocket from "znsocket";

import {
  FaArrowRotateRight,
  FaCircleInfo,
  FaFileCirclePlus,
  FaPencil,
} from "react-icons/fa6";
import { TbPlugConnected } from "react-icons/tb";
import { MdAddChart, MdExitToApp } from "react-icons/md";

function getServerUrl() {
  const { protocol, host } = window.location;
  const basePath = import.meta.env.BASE_URL || "/";
  return `${protocol}//${host}${basePath}`;
}

function ConsoleWindow({
  text,
  setConsoleShow,
  token,
  setConsoleText,
  colorMode,
  step,
  selection,
}: {
  text: string[];
  setConsoleShow: any;
  token: string;
  setConsoleText: any;
  colorMode: string;
  step: number;
  selection: Set<number>;
}) {
  const [showTime, setShowTime] = useState(false);
  const [chatInput, setChatInput] = useState<object>({});
  let [conInterface, setConInterface]: any = useState(undefined);
  const [showDropdown, setShowDropdown] = useState(false);
  let chatInputRef = useRef(null);



  const handleChatInputChange = (e: any) => {
    // log the selectionStart
    setChatInput({
      msg: e.target.value,
      time: new Date().toLocaleTimeString(),
    });
    if (e.target.value.endsWith('!!')) {
      setShowDropdown(true);
    } else {
      setShowDropdown(false);
    }
  };

  useEffect(() => {
    const con = new znsocket.List({
      client: client,
      key: "room:" + token + ":chat",
    });
    setConInterface(con);
  }, [token]);

  const handleSendMessage = () => {
    conInterface.append(chatInput);
    setConsoleText([...text, chatInput]);
    if (chatInputRef.current) {
      chatInputRef.current.value = "";
    }
  };

  const handleKeyPress = (e) => {
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      handleSendMessage();
    }
  };

  return (
    <>
    <Rnd
      default={{
        x: window.innerWidth / 2 - 400,
        y: window.innerHeight / 2 - 300,
        width: 380,
        height: 280,
      }}
      minHeight={200}
      minWidth={200}
      style={{
        zIndex: 1000,
        padding: 0,
        margin: 0,
      }}
    >
      <Card
        style={{
          margin: 0,
          padding: 0,
          width: "100%",
          height: "100%",
        }}
        // ref={cardRef}
      >
        <Card.Header className="d-flex justify-content-between align-items-center">
          <Card.Title>Console</Card.Title>
          <div className="d-flex align-items-center">
            <Form.Check
              type="switch"
              id="show-time-switch"
              label="Show Time"
              checked={showTime}
              onChange={() => {
                setShowTime(!showTime);
              }}
              className="me-2"
            />
            <Button variant="close" onClick={() => setConsoleShow(false)} />
          </div>
        </Card.Header>

        {/* Message Body with Optional Timestamp */}
        <Card.Body className="text-start overflow-y-auto">
          {text.map((line, idx) => (
            <p key={idx}>
              {showTime && <span className="text-muted me-2">{line.time}</span>}
              <Markdown
                remarkPlugins={[remarkMath, remarkGfm]}
                rehypePlugins={[rehypeKatex]}
                children={line.msg}
                components={{
                  code(props) {
                    const { children, className, node, ...rest } = props;
                    const match = /language-(\w+)/.exec(className || "");
                    return match ? (
                      <SyntaxHighlighter
                        {...rest}
                        PreTag="div"
                        children={String(children).replace(/\n$/, "")}
                        language={match[1]}
                        style={colorMode === "light" ? oneLight : oneDark}
                      />
                    ) : (
                      <code {...rest} className={className}>
                        {children}
                      </code>
                    );
                  },
                }}
              />
            </p>
          ))}
        </Card.Body>

        <Card.Footer>
          <InputGroup>
            <Form.Control
              as="textarea"
              rows={1}
              placeholder="Type a message..."
              // value={inputValue}
              // onChange={handleChatInputChange}
              onInput={handleChatInputChange}
              onKeyDown={handleKeyPress}
              ref={chatInputRef}
            />
            <Button variant="primary" onClick={handleSendMessage}>
              Send
            </Button>
          </InputGroup>
        </Card.Footer>
      </Card>
    </Rnd>
    {showDropdown && <ChatInsertModal show={showDropdown} onHide={() => setShowDropdown(false)} chatInputRef={chatInputRef} step={step} selection={selection}/>}
    </>
  );
}

function ChatInsertModal({ show, onHide, chatInputRef, step, selection }: any) {
  const options = [
    { value: 'step', label: 'step' },
    { value: 'selection', label: 'selection' },
  ];

  const handleSelectChange = (selectedOption: any) => {
    chatInputRef.current.value = chatInputRef.current.value.slice(0, -2)
    if (selectedOption.value === 'step') {
      chatInputRef.current.value = chatInputRef.current.value + `[${selectedOption.value}](${window.location.origin}/?step=${step})`;
    } else if (selectedOption.value === 'selection') {
      if (selection.size === 0) {
        chatInputRef.current.value = chatInputRef.current.value + `[${selectedOption.value}](${window.location.origin}/?selection=null)`;
      } else {
        chatInputRef.current.value = chatInputRef.current.value + `[${selectedOption.value}](${window.location.origin}/?selection=${Array.from(selection)})`;
      }
    }
    // trigger the change event
    chatInputRef.current.dispatchEvent(new InputEvent("input", { bubbles: true }));
    onHide();
  }

  return (
    <Modal show={show} aria-labelledby="contained-modal-title-vcenter" size="lg">
      <Modal.Header closeButton>
        <Modal.Title id="contained-modal-title-vcenter">
          ZnDraw Chat Insert
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>
        <Container>
          Share your current view:
        <Select
              options={options}
              onChange={handleSelectChange}
              // menuIsOpen={true}
              placeholder="Choose..."
            />
        </Container>
      </Modal.Body>
      <Modal.Footer>
        <Button onClick={onHide}>Close</Button>
      </Modal.Footer>
    </Modal>
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
- **add bookmark at current step**: \`keypress B\`
- **jump between bookmarks**: \`shift + keypress ▶\\◀\`
- **remove bookmark**: \`shift + mouse click\`
- **toggle drawing mode**: \`keypress X\`
- **select all particles**: \`ctrl + A\`
- **remove selected particles**: \`backspace\`
- **remove line point**: \`shift + backspace\`
- **reset camera to original position**: \`keypress O\`
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
  const serverUrl = getServerUrl();
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
          Python Client
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
  const serverUrl = getServerUrl();
  const urlWithRoom = `${serverUrl}token/${room}`;
  const resetURL = `${serverUrl}reset`;

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
        <Button href={resetURL}>Create new Scene</Button>
      </Modal.Footer>
    </Modal>
  );
}

interface TutorialModalProps {
  show: boolean;
  onHide: () => void;
  url: string;
}

const TutorialModal: React.FC<TutorialModalProps> = ({ show, onHide, url }) => {
  return (
    <Modal show={show} onHide={onHide} size="xl" dialogClassName="custom-modal">
      <Modal.Header closeButton>
        <Modal.Title>ZnDraw Tutorial</Modal.Title>
      </Modal.Header>
      <Modal.Body className="modal-body-custom">
        <iframe
          src={url}
          id="tutorialIframe"
          allowFullScreen
          className="iframe-custom"
        ></iframe>
      </Modal.Body>
    </Modal>
  );
};

function SiMGenButtons({ queuePosition }: { queuePosition: number }) {
  const runConnect = () => {
    socket.emit("modifier:run", {
      method: { discriminator: "Connect" },
    });
  };

  const runGenerate = () => {
    socket.emit("modifier:run", {
      method: { discriminator: "SiMGenDemo" },
    });
  };

  const createNewCanvas = () => {
    socket.emit("modifier:run", {
      method: { discriminator: "NewCanvas" },
    });
  };

  return (
    <>
      <BtnTooltip text="Connect selected atoms (shift click)">
        <Button
          variant="success"
          className="mx-1"
          onClick={runConnect}
          disabled={queuePosition != -1}
        >
          <TbPlugConnected /> Connect
        </Button>
      </BtnTooltip>
      <BtnTooltip text="Run SiMGen molecular generation">
        <Button
          variant="success"
          className="mx-1"
          onClick={runGenerate}
          disabled={queuePosition != -1}
        >
          <FaRocket /> Generate
        </Button>
      </BtnTooltip>
      <BtnTooltip text="Replace scene with empty canvas and enter drawing mode">
        <Button
          variant="success"
          className="mx-1"
          onClick={createNewCanvas}
          disabled={queuePosition != -1}
        >
          <FaFileCirclePlus /> New Canvas
        </Button>
      </BtnTooltip>
    </>
  );
}

const FileUpload = forwardRef((props, ref) => {
  const fileInputRef = useRef(null);

  const handleClick = () => {
    fileInputRef.current.click();
  };

  const handleFileChange = (event) => {
    const file = event.target.files[0];
    const formData = new FormData();
    formData.append("file", file);

    const basePath = import.meta.env.BASE_URL || "/";
    fetch(`${basePath}upload`, {
      method: "POST",
      body: formData,
    });

    if (ref) {
      // Check if ref is provided
      ref.current = fileInputRef.current; // Forward the ref to the underlying DOM element
    }
  };

  return (
    <div>
      <input
        type="file"
        ref={fileInputRef}
        onChange={handleFileChange}
        style={{ display: "none" }} // Hide the file input visually
      />
      <BtnTooltip text="Upload">
        <Button
          variant="outline-primary"
          className="mx-1"
          onClick={handleClick}
        >
          <FaUpload />
        </Button>
      </BtnTooltip>
    </div>
  );
});

interface HeadBarProps {
  room: string;
  colorMode: string;
  handleColorMode: any;
  setIsDrawing: any;
  setGeometries: any;
  setPoints: any;
  isDrawing: boolean;
  tutorialURL: string;
  showSiMGen: boolean;
  modifierQueue: number;
  isAuthenticated: boolean;
  roomLock: boolean;
  setAddPlotsWindow: any;
  messages: any[];
  token: string;
  step: number;
  selection: Set<number>;
}

const HeadBar = ({
  room,
  colorMode,
  handleColorMode,
  setIsDrawing,
  setGeometries,
  setPoints,
  isDrawing,
  tutorialURL,
  showSiMGen,
  modifierQueue,
  isAuthenticated,
  roomLock,
  setAddPlotsWindow,
  messages,
  token,
  step,
  selection,
}: HeadBarProps) => {
  const [helpModalShow, setHelpModalShow] = useState(false);
  const [connectModalShow, setConnectModalShow] = useState(false);
  const [refreshModalShow, setRefreshModalShow] = useState(false);
  const [tutorialModalShow, setTutorialModalShow] = useState(false);
  const [consoleShow, setConsoleShow] = useState(false);
  const [consoleText, setConsoleText] = useState<string[]>([]);

  useEffect(() => {
    setConsoleText(messages);
  }, [messages]);

  useEffect(() => {
    setConsoleShow(showSiMGen);
  }, [showSiMGen]);

  const handleRemovePointsGeometries = () => {
    console.log("remove points and geometries");
    socket.emit("room:geometry:set", []);
    socket.emit("room:point:set", []);
    setGeometries([]);
    setPoints([]);
  };
  const basePath = useMemo(() => import.meta.env.BASE_URL || "/", []);

  return (
    <>
      <Navbar
        expand="md"
        className="bg-body-tertiary fixed-top"
        style={{ height: 50 }}
      >
        <Container fluid>
          <Navbar.Brand>
            <Button
              className="px-0 py-0 btn-lg"
              variant="tertiary"
              href="https://github.com/zincware/zndraw"
              target="_blank"
            >
              ZnDraw
            </Button>
            {showSiMGen && (
              <>
                <Button className="px-0 py-0 btn-lg" variant="tertiary">
                  +
                </Button>
                <Button
                  className="px-0 py-0 btn-lg"
                  variant="tertiary"
                  href="https://github.com/zincware/zndraw"
                  target="_blank"
                >
                  SiMGen
                </Button>
              </>
            )}
          </Navbar.Brand>
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
              {showSiMGen && <SiMGenButtons queuePosition={modifierQueue} />}
            </Nav>
            <Nav className="ms-auto">
              {tutorialURL && (
                <BtnTooltip text="Tutorial">
                  <Button
                    variant="warning"
                    className="mx-1"
                    onClick={() => setTutorialModalShow(true)}
                  >
                    Tutorial <FaFilm />
                  </Button>
                </BtnTooltip>
              )}
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
              <BtnTooltip text="Add plots window">
                <Button
                  variant="outline-primary"
                  className="mx-1"
                  onClick={() => setAddPlotsWindow((prev: number) => prev + 1)}
                >
                  <MdAddChart />
                </Button>
              </BtnTooltip>
              <FileUpload />
              <BtnTooltip text="Download">
                <Button
                  variant="outline-primary"
                  className="mx-1"
                  href={`${basePath}download`}
                  target="_blank"
                >
                  <FaDownload />
                </Button>
              </BtnTooltip>
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
              {isAuthenticated && (
                <>
                  <BtnTooltip
                    text={roomLock ? "Unlock this room" : "Lock this room"}
                  >
                    <Button
                      variant="outline-danger"
                      className="mx-1"
                      onClick={() => {
                        socket.emit("room:lock:set", !roomLock);
                      }}
                    >
                      {roomLock ? <FaLock /> : <FaLockOpen />}
                    </Button>
                  </BtnTooltip>
                  {/* <Button variant="outline-primary" className="mx-1">
                <FaUsers />
              </Button> */}
                  <BtnTooltip text="Close ZnDraw">
                    <Button
                      variant="outline-danger"
                      className="mx-1"
                      onClick={() => {
                        socket.emit("shutdown");
                        close();
                      }}
                    >
                      <MdExitToApp />
                    </Button>
                  </BtnTooltip>
                </>
              )}
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
      <TutorialModal
        show={tutorialModalShow}
        onHide={() => setTutorialModalShow(false)}
        url={tutorialURL}
      />
      {consoleShow && (
        <ConsoleWindow
          text={consoleText}
          setConsoleShow={setConsoleShow}
          token={token}
          setConsoleText={setConsoleText}
          colorMode={colorMode}
          step={step}
          selection={selection}
        />
      )}
    </>
  );
};

export default HeadBar;
