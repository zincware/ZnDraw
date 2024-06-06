import React from "react";
import {
  ProgressBar,
  InputGroup,
  Form,
  Container,
  Row,
  Col,
} from "react-bootstrap";

import { FaRegBookmark } from "react-icons/fa";
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
              placement="top"
              delay={{ show: 0, hide: 100 }}
              overlay={<Tooltip>{name}</Tooltip>}
            >
              <Col
                xs="auto"
                key={number}
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
  bookmarks: { number: string };
  setBookmarks: any;
}

const FrameProgressBar: React.FC<FrameProgressBarProps> = ({
  length,
  step,
  setStep,
  bookmarks,
  setBookmarks,
}) => {
  const handleClick = (event: React.MouseEvent<HTMLDivElement, MouseEvent>) => {
    const percentage = event.clientX / window.innerWidth;
    const newStep = Math.floor(percentage * length);
    setStep(newStep);
  };

  return (
    <>
      <Container
        fluid
        className="fixed-bottom py-0"
        style={{ pointerEvents: "none" }}
      >
        <Row
          className="justify-content-center"
          style={{ pointerEvents: "none" }}
        >
          <Col xs="auto" style={{ pointerEvents: "auto" }}>
            <JumpFrame step={step} setStep={setStep} length={length} />
          </Col>
        </Row>
        <Row
          className="justify-content-center"
          style={{ pointerEvents: "auto" }}
        >
          <Col xs={12} className="px-0 py-0">
            <ProgressBar
              now={step}
              className="frame-progress-bar"
              onClick={handleClick}
              max={length - 1}
              style={{ margin: "0 auto", padding: 0, cursor: "pointer" }}
            />
          </Col>
        </Row>
      </Container>
      <Bookmarks
        bookmarks={bookmarks}
        length={length}
        setStep={setStep}
        setBookmarks={setBookmarks}
      />
    </>
  );
};

export default FrameProgressBar;
