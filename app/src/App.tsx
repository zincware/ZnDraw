import React from "react";
import { AppProvider } from "./contexts/AppContext";
import { VisualizationContainer } from "./components/containers/VisualizationContainer";
import { UIContainer } from "./components/containers/UIContainer";
import { useSocketManager } from "./hooks/useSocketManager";
import { useKeyboardHandler } from "./hooks/useKeyboardHandler";
import { useFileHandler } from "./hooks/useFileHandler";
import { useSelectionCleanup } from "./hooks/useSelectionCleanup";
import { useColorMode } from "./components/utils";
import "./App.css";

// Internal component that uses all the hooks
const AppContent: React.FC = () => {
	// Initialize all hooks
	useSocketManager();
	useKeyboardHandler();
	// useWheelHandler();
	useSelectionCleanup();

	const { onDragOver, onDrop, onPointerMissed } = useFileHandler();

	return (
		<>
			<VisualizationContainer
				onPointerMissed={onPointerMissed}
				onDragOver={onDragOver}
				onDrop={onDrop}
			/>
			<UIContainer />
		</>
	);
};

export default function App() {
	const [colorMode, handleColorMode] = useColorMode();

	return (
		<AppProvider colorMode={colorMode} handleColorMode={handleColorMode}>
			<AppContent />
		</AppProvider>
	);
}
