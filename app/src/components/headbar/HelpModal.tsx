import {
	Button,
	Container,
	Dialog,
	DialogActions,
	DialogContent,
	DialogTitle,
} from "@mui/material";
import type React from "react";
import Markdown from "react-markdown";

interface HelpModalProps {
	show: boolean;
	onHide: () => void;
}

export function HelpModel(props: HelpModalProps) {
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
- **rotate camera perpendicular to the screen**: \`keypress R\` or \`keypress ctrl + R\`
- (experimental) **Enter molecule editor (none, translate, rotate)**: \`keypress E\`
`;

	return (
		<Dialog
			{...props}
			aria-labelledby="contained-modal-title-vcenter"
			maxWidth="lg"
			fullWidth
			open={props.show}
			onClose={props.onHide}
		>
			<DialogTitle id="contained-modal-title-vcenter">ZnDraw Help</DialogTitle>

			<DialogContent>
				<Container>
					<Markdown>{helpMD}</Markdown>
				</Container>
			</DialogContent>
			<DialogActions>
				<Button onClick={props.onHide}>Close</Button>
			</DialogActions>
		</Dialog>
	);
}
