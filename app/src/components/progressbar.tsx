import React, { useState, useEffect } from "react";
import { InputGroup, Form, Container, Row, Col, Card } from "react-bootstrap";

import { FaEye, FaLock, FaRegBookmark } from "react-icons/fa";
import { OverlayTrigger, Tooltip } from "react-bootstrap";
import { PiSelectionSlash } from "react-icons/pi";

interface JumpFrameProps {
  step: number;
  setStep: (step: number) => void;
  length: number;
}

const JumpFrame: React.FC<JumpFrameProps> = ({ step, setStep, length }) => {
  const handleBlur = (e: React.FocusEvent<HTMLInputElement>) => {
    if (e.target.value === "") {
      return;
    }
    const newStep = parseInt(e.target.value, 10);
    if (newStep >= 0 && newStep < length) {
      setStep(newStep);
    } else {
      alert(
        "Invalid input. Please enter a number between 0 and " + (length - 1),
      );
    }
    e.target.value = "";
  };

  const handleKeyDown = (e: React.KeyboardEvent<HTMLInputElement>) => {
    if (e.key === "Enter") {
      e.preventDefault();
      e.currentTarget.blur();
    }
  };

  return (
    <InputGroup>
      <Form.Control
        className="text-center"
        placeholder={`${step}/${length - 1}`}
        onBlur={handleBlur}
        onKeyDown={handleKeyDown}
        style={{
          background: "transparent",
          borderColor: "transparent",
          zIndex: 1,
        }}
      />
    </InputGroup>
  );
};

interface ProgressBarProps {
  length: number;
  disabledFrames: number[];
  bookmarks: any;
  setBookmarks: any;
  step: number;
  setStep: (step: number) => void;
}

const ColoredTiles = ({
  length,
  disabledFrames,
  setStep,
  tickInterval,
}: {
  length: number;
  disabledFrames: number[];
  setStep: (step: number) => void;
  tickInterval: number;
}) => {
  useEffect(() => {
    console.log("component rerendered");
  }, [length, disabledFrames, setStep, tickInterval]);

  const [disabledPositions, setdisabledPositions] = useState<number[]>([]);
  const [ticks, setTicks] = useState<number[]>([]);

  useEffect(() => {
    const disabledPositions = [...Array(length).keys()].filter((position) =>
      disabledFrames.includes(position),
    );
    setdisabledPositions(disabledPositions);
  }, [length, disabledFrames]);

  useEffect(() => {
    const ticks = [...Array(length).keys()].filter(
      (position) => position % tickInterval === 0,
    );
    setTicks(ticks);
  }, [length, tickInterval]);

  const onTileClick = (event: any) => {
    const rect = event.target.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const position = Math.floor((x / rect.width) * length);
    setStep(position);
  };
  return (
    <>
      {ticks.map((position) => {
        const commonStyles = {
          left: `${(position / length) * 100}%`,
          height: 25,
        };
        return (
          <div
            key={position}
            className={`position-absolute`}
            style={commonStyles}
          >
            <div className="progress-bar-tick-line bg-dark"></div>
          </div>
        );
      })}
      <div
        className={`position-absolute bg-gradient bg-primary`}
        style={{ width: "100%", height: 25 }}
        onClick={(e) => onTileClick(e)}
      ></div>

      {disabledPositions.map((position) => {
        const commonStyles = {
          left: `${(position / length) * 100}%`,
          width: `${100 / (length - 1)}%`,
          height: 25,
        };
        return (
          <div
            key={position}
            className={`position-absolute p-0 bg-gradient bg-primary-subtle`}
            style={commonStyles}
          ></div>
        );
      })}
    </>
  );
};

const Bookmarks = ({
  length,
  bookmarks,
  setBookmarks,
  setStep,
}: {
  length: number;
  bookmarks: { [key: number]: string };
  setBookmarks: (bookmarks: { [key: number]: string }) => void;
  setStep: (step: number) => void;
}) => {
  const handleBookmarkClick = (event: any, number: number) => {
    if (event.shiftKey) {
      const newBookmarks = { ...bookmarks };
      delete newBookmarks[number];
      setBookmarks(newBookmarks);
    } else {
      setStep(number);
    }
  };

  return (
    <>
      {Object.keys(bookmarks).map((key) => {
        const position = parseInt(key);
        return (
          <OverlayTrigger
            key={position}
            placement="top"
            delay={{ show: 0, hide: 100 }}
            overlay={
              <Tooltip style={{ marginLeft: "0.375em" }}>
                {bookmarks[position]}
              </Tooltip>
            }
          >
            <div
              className="position-absolute progress-bar-bookmark"
              style={{
                left: `${(position / length) * 100}%`,
                top: 7,
                marginLeft: "-0.375em",
              }}
            >
              <FaRegBookmark
                className="position-absolute"
                size={"0.75em"}
                onClick={(e) => handleBookmarkClick(e, position)}
              />
            </div>
          </OverlayTrigger>
        );
      })}
    </>
  );
};

const VLine = ({ length, step }: { length: number; step: number }) => {
  return (
    <div
      className="position-absolute p-0 progress-bar-v-line"
      style={{ left: `${(step / length) * 100}%`, height: 25, width: "2px" }}
      // why do I need to overwrite the height and width here?
    ></div>
  );
};

const ProgressBar = ({
  length,
  disabledFrames,
  bookmarks,
  setBookmarks,
  step,
  setStep,
}: ProgressBarProps) => {
  const [tickInterval, setTickInterval] = useState<number>(1);

  useEffect(() => {
    setTickInterval(Math.floor(length / 100) + 1);
  }, [length]);

  return (
    <Row className="position-relative">
      <ColoredTiles
        length={length}
        disabledFrames={disabledFrames}
        setStep={setStep}
        tickInterval={tickInterval}
      />
      <Bookmarks
        length={length}
        bookmarks={bookmarks}
        setBookmarks={setBookmarks}
        setStep={setStep}
      />
      <VLine length={length} step={step} />
    </Row>
  );
};

interface FrameProgressBarProps {
  length: number;
  step: number;
  setStep: (step: number) => void;
  selectedFrames: Set<number>;
  bookmarks: any[]; // Replace with actual type if known
  setBookmarks: (bookmarks: any[]) => void; // Replace with actual type if known
  setSelectedFrames: (selectedFrames: Set<number>) => void;
}

const FrameProgressBar: React.FC<FrameProgressBarProps> = ({
  length,
  step,
  setStep,
  selectedFrames,
  bookmarks,
  setBookmarks,
  setSelectedFrames,
}) => {
  const [linePosition, setLinePosition] = useState<number>(0);
  const [disabledFrames, setDisabledFrames] = useState<number[]>([]);
  const progressHandleParentRef = React.createRef<HTMLDivElement>();

  useEffect(() => {
    // disable frames are the frames that are not selected, if the selectedFrames is not empty
    if (selectedFrames.size > 0) {
      const disabledFrames = [...Array(length).keys()].filter(
        (frame) => !selectedFrames.has(frame),
      );
      setDisabledFrames(disabledFrames);
    } else {
      setDisabledFrames([]);
    }
  }, [selectedFrames, length]);

  const handleSelectionReset = () => {
    console.log("Resetting selection");
    setSelectedFrames(new Set());
  };

  useEffect(() => {
    // Calculate the linePosition based on the step, length, and window width
    setLinePosition((step / length) * 100);
  }, [step, length]);

  const handleMouseDown = (e) => {
    e.preventDefault();
    if (!progressHandleParentRef.current) {
      console.error(
        "progressHandleParentRef is not set - dragging will not work",
      );
      return;
    }
    // we need to store the parent rect in a variable because
    // once we move the mouse, the parent rect will change
    const parentRect = progressHandleParentRef.current.getBoundingClientRect();

    const handleMouseMove = (e) => {
      // compute the lineposition based on the mouse position and the length
      const x = e.clientX - parentRect.left;
      const position = Math.floor((x / parentRect.width) * length);
      setStep(Math.max(0, position));
    };

    document.addEventListener("mousemove", handleMouseMove);

    document.addEventListener(
      "mouseup",
      () => {
        document.removeEventListener("mousemove", handleMouseMove);
      },
      { once: true },
    );
  };

  return (
    <Container fluid className="fixed-bottom px-0 py-0">
      <Row>
        <Col xs="2">
          <Row>
            <Col
              className="d-flex bg-dark bg-gradient justify-content-center align-items-center"
              style={{ height: 1 }}
            ></Col>
          </Row>
          <Row>
            <Col
              className="d-flex bg-secondary-subtle bg-gradient justify-content-center align-items-center"
              style={{ height: 25 }}
            >
              <JumpFrame step={step} setStep={setStep} length={length} />
              <OverlayTrigger
                placement="top"
                delay={{ show: 0, hide: 100 }}
                overlay={<Tooltip>reset selection</Tooltip>}
              >
                <div>
                  {" "}
                  <PiSelectionSlash onClick={handleSelectionReset} />
                </div>
              </OverlayTrigger>
            </Col>
          </Row>
        </Col>
        <Col>
          <Row className="position-relative">
            <Col
              className="d-flex justify-content-center"
              ref={progressHandleParentRef}
            >
              <div
                className="handle"
                style={{ left: `${linePosition}%`, cursor: "pointer" }}
                onMouseDown={(e) => handleMouseDown(e)}
              >
                <div className="square"></div>
                <div className="triangle"></div>
              </div>
            </Col>
          </Row>
          <Row>
            <Col
              className="d-flex bg-dark bg-gradient justify-content-center align-items-center"
              style={{ height: 1 }}
            ></Col>
          </Row>

          <ProgressBar
            length={length}
            disabledFrames={disabledFrames}
            bookmarks={bookmarks}
            step={step}
            setStep={setStep}
            setBookmarks={setBookmarks}
          />
        </Col>
      </Row>
    </Container>
  );
};

export default FrameProgressBar;
