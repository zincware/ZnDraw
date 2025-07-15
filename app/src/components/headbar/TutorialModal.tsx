import { Dialog, DialogContent, DialogTitle } from "@mui/material";
import type React from "react";

interface TutorialModalProps {
	show: boolean;
	onHide: () => void;
	url: string;
}

export const TutorialModal: React.FC<TutorialModalProps> = ({
	show,
	onHide,
	url,
}) => {
	return (
		<Dialog open={show} onClose={onHide} maxWidth="xl" fullWidth>
			<DialogTitle>ZnDraw Tutorial</DialogTitle>
			<DialogContent>
				<iframe
					src={url}
					id="tutorialIframe"
					allowFullScreen
					style={{ width: "100%", height: "80vh", border: "none" }}
				/>
			</DialogContent>
		</Dialog>
	);
};
