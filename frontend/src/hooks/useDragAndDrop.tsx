import { type DragEvent, useCallback, useState } from "react";

interface UseDragAndDropReturn {
	isDragging: boolean;
	handleDragOver: (e: DragEvent<HTMLDivElement>) => void;
	handleDragEnter: (e: DragEvent<HTMLDivElement>) => void;
	handleDragLeave: (e: DragEvent<HTMLDivElement>) => void;
	handleDrop: (e: DragEvent<HTMLDivElement>) => void;
}

export function useDragAndDrop(
	onFiles: (files: File[]) => Promise<void>,
): UseDragAndDropReturn {
	const [isDragging, setIsDragging] = useState(false);
	const [dragCounter, setDragCounter] = useState(0);

	const handleDragOver = useCallback((e: DragEvent<HTMLDivElement>) => {
		e.preventDefault();
		e.stopPropagation();
	}, []);

	const handleDragEnter = useCallback((e: DragEvent<HTMLDivElement>) => {
		e.preventDefault();
		e.stopPropagation();
		setDragCounter((prev) => prev + 1);
		setIsDragging(true);
	}, []);

	const handleDragLeave = useCallback((e: DragEvent<HTMLDivElement>) => {
		e.preventDefault();
		e.stopPropagation();
		setDragCounter((prev) => {
			const newCount = prev - 1;
			if (newCount === 0) {
				setIsDragging(false);
			}
			return newCount;
		});
	}, []);

	const handleDrop = useCallback(
		async (e: DragEvent<HTMLDivElement>) => {
			e.preventDefault();
			e.stopPropagation();
			setIsDragging(false);
			setDragCounter(0);

			if (!e.dataTransfer.files || e.dataTransfer.files.length === 0) {
				return;
			}

			await onFiles(Array.from(e.dataTransfer.files));
		},
		[onFiles],
	);

	return {
		isDragging,
		handleDragOver,
		handleDragEnter,
		handleDragLeave,
		handleDrop,
	};
}
