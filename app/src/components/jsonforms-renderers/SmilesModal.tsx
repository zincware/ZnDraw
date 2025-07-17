import {
	useState,
} from "react";
import {
	Dialog,
	DialogContent,
	DialogActions,
	Button,
	Box,
} from "@mui/material";
import { Jsme } from "@loschmidt/jsme-react";

export const SmilesModal = ({
	smiles,
	setSmiles,
	onClose, // Recommended: allow parent to control dialog visibility
}: {
	smiles: string;
	setSmiles: (newSmiles: string) => void;
	onClose: () => void;
}) => {
	const [editorSmiles, setEditorSmiles] = useState(smiles || "");

	const handleSave = () => {
		if (editorSmiles) {
			setSmiles(editorSmiles);
			console.log("Saved SMILES:", editorSmiles);
			onClose?.(); // Close modal after saving
		}
	};

	return (
		<Dialog open={true} onClose={onClose} fullWidth maxWidth="sm">
			<DialogContent dividers>
				<Box display="flex" justifyContent="center">
					<Jsme
						height={500}
						width={500}
						options=""
						onChange={setEditorSmiles}
						src="https://jsme-editor.github.io/dist/jsme/jsme.nocache.js"
						smiles={smiles}
					/>
				</Box>
			</DialogContent>

			<DialogActions>
				<Button onClick={onClose}>Cancel</Button>
				<Button
					onClick={handleSave}
					variant="contained"
					color="primary"
					disabled={!editorSmiles}
				>
					Save
				</Button>
			</DialogActions>
		</Dialog>
	);
};

export default SmilesModal;