import { Box } from "@mui/material";
import { useColorScheme } from "@mui/material/styles";
import type React from "react";
import { UIContainer } from "./components/containers/UIContainer";
import { VisualizationContainer } from "./components/containers/VisualizationContainer";
import { AppProvider } from "./contexts/AppContext";
import { SlowFrameProvider } from "./contexts/SlowFrameContext";
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
			<UIContainer />
		</Box>
	);
};

export default function App() {
	const { mode, setMode } = useColorScheme();

	const handleColorMode = () => {
		const newMode = mode === "light" ? "dark" : "light";
		setMode(newMode);
	};

	return (
		<AppProvider colorMode={mode || "light"} handleColorMode={handleColorMode}>
			<SlowFrameProvider>
				<AppContent />
			</SlowFrameProvider>
		</AppProvider>
	);
}
