import { Card } from "react-bootstrap";

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
