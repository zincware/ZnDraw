import type React from "react";
import { UIContainer } from "./components/containers/UIContainer";
import { VisualizationContainer } from "./components/containers/VisualizationContainer";
import { useColorMode } from "./components/utils";
import { AppProvider } from "./contexts/AppContext";
import { useFileHandler } from "./hooks/useFileHandler";
import { useKeyboardHandler } from "./hooks/useKeyboardHandler";
import { useSelectionCleanup } from "./hooks/useSelectionCleanup";
import { useSocketManager } from "./hooks/useSocketManager";
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
