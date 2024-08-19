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

const ProgressBar = ({
  length,
  disabledFrames,
  bookmarks,
  setBookmarks,
  step,
  setStep,
}: ProgressBarProps) => {
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
    <Row className="position-relative">
      {[...Array(length).keys()].map((position) => {
        const isDisabled = disabledFrames.includes(position);
        const commonStyles = {
          left: `${(position / length) * 100}%`,
          height: 25,
        };

        return (
          <div
            key={position}
            className={`position-absolute bg-gradient ${isDisabled ? "bg-primary-subtle" : "bg-primary"}`}
            style={commonStyles}
            onClick={isDisabled ? undefined : () => setStep(position)}
          >
            {bookmarks[position] && (
              <OverlayTrigger
                placement="top"
                delay={{ show: 0, hide: 100 }}
                overlay={<Tooltip>{bookmarks[position]}</Tooltip>}
              >
                {/* overlay trigger requires the div wrapped around the icon */}
                <div className="progress-bar-bookmark">
                  <FaRegBookmark
                    size={"0.75em"}
                    className="progress-bar-bookmark"
                    onClick={(e) => handleBookmarkClick(e, position)}
                  />
                </div>
              </OverlayTrigger>
            )}
            <div className="progress-bar-tick-line bg-dark"></div>
            {position === step && <div className="progress-bar-v-line"></div>}
          </div>
        );
      })}
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

  return (
    <Container fluid className="fixed-bottom px-0 py-0">
      <Row>
        <Col xs="1">
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
              <OverlayTrigger placement="top" delay={{ show: 0, hide: 100 }} overlay={<Tooltip>reset selection</Tooltip>}>
              <div>                <PiSelectionSlash onClick={handleSelectionReset} />
              </div>
              </OverlayTrigger>
            </Col>
          </Row>
        </Col>
        <Col>
          <Row className="position-relative">
            <Col className="d-flex justify-content-center">
              <div className="handle" style={{ left: `${linePosition}%` }}>
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
