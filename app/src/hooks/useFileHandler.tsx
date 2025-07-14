import { useCallback } from "react";
import { socket } from "../socket";
import { useAppContext } from "../contexts/AppContext";

export const useFileHandler = () => {
	const { setSelectedPoint, setSelectedIds, isDrawing, setIsDrawing } =
		useAppContext();

	const onDragOver = useCallback((event: React.DragEvent) => {
		event.preventDefault();
	}, []);

	const onDrop = useCallback(async (event: React.DragEvent) => {
		event.preventDefault();

		const file = event.dataTransfer.files[0];
		if (!file) {
			console.error("No file was dropped");
			return;
		}

		try {
			const arrayBuffer = await file.arrayBuffer();
			const content = new Uint8Array(arrayBuffer);

			// send the file to the server
			socket.emit("room:upload:file", {
				content: Array.from(content),
				filename: file.name,
			});
		} catch (error) {
			console.error("Error reading file:", error);
		}
	}, []);

	const onPointerMissed = useCallback(() => {
		setSelectedPoint(null);
		setSelectedIds(new Set());
		if (isDrawing) {
			setIsDrawing(false);
		}
	}, [setSelectedPoint, setSelectedIds, isDrawing, setIsDrawing]);

	return {
		onDragOver,
		onDrop,
		onPointerMissed,
	};
};
