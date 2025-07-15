import { Button } from "@mui/material";
import { forwardRef, useRef } from "react";
import { FaUpload } from "react-icons/fa";
import { BtnTooltip } from "../tooltips";

export const FileUpload = forwardRef((props, ref) => {
	const fileInputRef = useRef(null);

	const handleClick = () => {
		fileInputRef.current.click();
	};

	const handleFileChange = (event) => {
		const file = event.target.files[0];
		const formData = new FormData();
		formData.append("file", file);

		const basePath = import.meta.env.BASE_URL || "/";
		fetch(`${basePath}upload`, {
			method: "POST",
			body: formData,
		});

		if (ref) {
			// Check if ref is provided
			ref.current = fileInputRef.current; // Forward the ref to the underlying DOM element
		}
	};

	return (
		<div>
			<input
				type="file"
				ref={fileInputRef}
				onChange={handleFileChange}
				style={{ display: "none" }} // Hide the file input visually
			/>
			<BtnTooltip text="Upload">
				<Button
					variant="outline-primary"
					className="mx-1"
					onClick={handleClick}
				>
					<FaUpload />
				</Button>
			</BtnTooltip>
		</div>
	);
});
