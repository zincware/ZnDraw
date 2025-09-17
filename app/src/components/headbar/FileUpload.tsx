import { forwardRef, useRef } from "react";
import type { ReactElement } from "react";

// The component now takes a function to render the button.
interface FileUploadProps {
	renderButton: (handleClick: () => void) => ReactElement;
}

export const FileUpload = forwardRef<HTMLInputElement, FileUploadProps>(
	({ renderButton }, ref) => {
		const fileInputRef = useRef<HTMLInputElement>(null);

		// This function will be passed to the renderButton prop.
		const handleClick = () => {
			fileInputRef.current?.click();
		};

		const handleFileChange = (event: React.ChangeEvent<HTMLInputElement>) => {
			const file = event.target.files?.[0];
			if (!file) return;

			const formData = new FormData();
			formData.append("file", file);

			const basePath = import.meta.env.BASE_URL || "/";
			fetch(`${basePath}upload`, {
				method: "POST",
				body: formData,
			});

			if (ref && typeof ref !== "function") {
				ref.current = fileInputRef.current;
			}
		};

		return (
			<>
				{/* The actual file input remains hidden */}
				<input
					type="file"
					ref={fileInputRef}
					onChange={handleFileChange}
					style={{ display: "none" }}
				/>
				{/* Call the render function, passing it the click handler */}
				{renderButton(handleClick)}
			</>
		);
	},
);

FileUpload.displayName = "FileUpload";
