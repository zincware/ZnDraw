import { Button, Card } from "react-bootstrap";
import { Rnd } from "react-rnd";
import type { Frame } from "./particles";

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
								<>
									<strong>{key}: </strong> {value}
									<br />
								</>
							))}
						</Card.Text>
					</Card.Body>
				</Card>
			)}
		</>
	);
};

export const SceneInfoOverlay = ({
	frame,
	setShowParticleInfo,
}: {
	frame: Frame;
	setShowParticleInfo: any;
}) => {
	return (
		<Rnd
			default={{
				x: window.innerWidth / 2 - 300,
				y: -window.innerHeight / 2 + 75,
				width: 280,
				height: "100px",
			}}
			minHeight={150}
			minWidth={150}
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
			>
				<Card.Header className="d-flex justify-content-between align-items-center">
					<Card.Title>Info</Card.Title>
					<Button variant="close" onClick={() => setShowParticleInfo(false)} />
				</Card.Header>
				<Card.Body className="text-start overflow-y-auto">
					{frame.calc.energy && (
						<>
							Energy: {frame.calc.energy} eV
							<br />
						</>
					)}
					Particles: {frame.positions.length}
				</Card.Body>
			</Card>
		</Rnd>
	);
};
