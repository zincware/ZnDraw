import React, { useState, useEffect } from "react";
import {
  ProgressBar,
  InputGroup,
  Form,
  Container,
  Row,
  Col,
  Card,
} from "react-bootstrap";

import { FaEye, FaLock, FaRegBookmark } from "react-icons/fa";
import { OverlayTrigger, Tooltip } from "react-bootstrap";

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
    <InputGroup className="mb-1">
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

const Bookmarks = ({
  setBookmarks,
  bookmarks,
  length,
  setStep,
}: {
  setBookmarks: any;
  bookmarks: { [number: string]: string };
  length: number;
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
    <Container fluid className="fixed-bottom py-0">
      <Row className="justify-content-center position-relative">
        {Object.entries(bookmarks).map(([number, name]) => {
          const position = (parseInt(number) / length) * 100;
          return (
            <OverlayTrigger
              key={`bookmark-${number}`}
              placement="top"
              delay={{ show: 0, hide: 100 }}
              overlay={<Tooltip>{name}</Tooltip>}
            >
              <Col
                xs="auto"
                style={{
                  position: "absolute",
                  left: `${position}%`,
                  bottom: "14px",
                  transform: "translateX(-50%)",
                }}
              >
                <FaRegBookmark
                  title={name}
                  color="grey"
                  onClick={(e) => handleBookmarkClick(e, parseInt(number))}
                />
              </Col>
            </OverlayTrigger>
          );
        })}
      </Row>
    </Container>
  );
};

interface FrameProgressBarProps {
  length: number;
  step: number;
  setStep: (step: number) => void;
  selectedFrames: number[];
  bookmarks: any[]; // Replace with actual type if known
  setBookmarks: (bookmarks: any[]) => void; // Replace with actual type if known
}

const FrameProgressBar: React.FC<FrameProgressBarProps> = ({
  length,
  step,
  setStep,
  selectedFrames,
  bookmarks,
  setBookmarks,
}) => {
  const [linePosition, setLinePosition] = useState<number>(0);
  const [disabledFrames, setDisabledFrames] = useState<number[]>([]);

  useEffect(() => {
    // disable frames are the frames that are not selected, if the selectedFrames is not empty
    if (selectedFrames.length > 0) {
      const disabledFrames = [...Array(length).keys()].filter(
        (frame) => !selectedFrames.includes(frame),
      );
      setDisabledFrames(disabledFrames);
    }
  }, [selectedFrames, length]);

  useEffect(() => {
    // Calculate the linePosition based on the step, length, and window width
    setLinePosition((step / (length)) * 100);
  }, [step, length]);

  return (
    <Container fluid className="fixed-bottom px-0 py-0">
      <Row>
        <Col xs="1">
          <Row>
            <Col
              className="d-flex bg-info justify-content-center align-items-center"
              style={{ height: 25 }}
            >
              <FaLock />
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

          <Row className="position-relative">
            {/* we add the last frame to click on */}
            {[...Array(length).keys()].map((position) =>
              disabledFrames.includes(position) ? (
                <div
                  key={position}
                  className="bg-secondary position-absolute" // A different class for disabled frames
                  style={{
                    left: `${(position / (length)) * 100}%`,
                    height: 25,
                  }}
                />
              ) : (
                <div
                  key={position}
                  className="bg-primary position-absolute"
                  style={{
                    left: `${(position / (length)) * 100}%`,
                    height: 25,
                  }}
                  onClick={() => setStep(position)}
                />
              ),
            )}
          </Row>
        </Col>
        {/* <Col xs="1">
          <Row>
            <Col
              className="d-flex bg-dark justify-content-center align-items-center"
              style={{ height: 25 }}
            >
              <FaEye />
            </Col>
          </Row>
        </Col> */}
      </Row>
    </Container>
  );
};

export default FrameProgressBar;
