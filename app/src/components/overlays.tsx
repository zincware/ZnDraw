import { Card } from "react-bootstrap";
import { Frame } from "./particles";
import { Rnd } from "react-rnd";

export const ParticleInfoOverlay = ({
  show,
  info,
  position,
}: {
  show: boolean;
  info: { [key: string]: any };
  position: { x: number; y: number };
}) => {
  return (
    <>
      {show && (
        <Card
          style={{
            position: "absolute",
            top: position.y + 5,
            left: position.x + 5,
            zIndex: 1000,
            padding: 0,
            margin: 0,
            maxWidth: "18rem",
          }}
        >
          <Card.Body>
            <Card.Text className="text-start">
              {Object.entries(info).map(([key, value]) => (
                <p key={key}>
                  <strong>{key}: </strong> {value}
                </p>
              ))}
            </Card.Text>
          </Card.Body>
        </Card>
      )}
    </>
  );
};

export const SceneInfoOverlay = ({ frame }: { frame: Frame }) => {
  console.log(frame);
  return (
    <Rnd
      default={{
        x: window.innerWidth / 2 - 300,
        y: -window.innerHeight / 2 + 75,
        width: 280,
        height: "100px",
      }}
      style={{ zIndex: 1000, padding: 0, margin: 0 }}
      i
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
          <Card.Title>Info</Card.Title>
        </Card.Header>
        <Card.Body>
          <Card.Text className="text-start text-nowrap">
            {frame.calc["energy"] && (
              <>
                Energy: {frame.calc["energy"]} eV
                <br />
              </>
            )}
            Particles: {frame.positions.length}
          </Card.Text>
        </Card.Body>
      </Card>
    </Rnd>
  );
};
