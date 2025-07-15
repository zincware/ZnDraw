import { Box } from "@mui/material";
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
		<Box sx={{ position: "relative", width: "100vw", height: "100vh" }}>
			<VisualizationContainer
				onPointerMissed={onPointerMissed}
				onDragOver={onDragOver}
				onDrop={onDrop}
			/>
			<Box
				s={{ position: "absolute", top: 0, left: 0, width: "100%", height: "100%", pointerEvents: "none" }}
			>
				<UIContainer />
			</Box>
		</Box>
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
