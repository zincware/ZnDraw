import { useEffect } from "react";
import { useAppContext } from "../contexts/AppContext";

export const useSelectionCleanup = () => {
	const { currentFrame, selectedIds, setSelectedIds } = useAppContext();

	useEffect(() => {
		if (currentFrame.positions.length > 0 && selectedIds.size > 0) {
			const newSelectedIds = new Set(
				Array.from(selectedIds).filter(
					(id) => id < currentFrame.positions.length,
				),
			);
			// if the selection is reduced, update the selection
			if (newSelectedIds.size < selectedIds.size) {
				setSelectedIds(newSelectedIds);
			}
		}
	}, [currentFrame, selectedIds, setSelectedIds]);
};
