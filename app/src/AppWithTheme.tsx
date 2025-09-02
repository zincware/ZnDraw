import { CssBaseline } from "@mui/material";
import {
	ThemeProvider,
	createTheme,
	useColorScheme,
} from "@mui/material/styles";
import type React from "react";
import { useEffect, useState } from "react";
import App from "./App";

const theme = createTheme({
	cssVariables: {
		colorSchemeSelector: "data-mui-color-scheme",
	},
	defaultColorScheme: "light",
	colorSchemes: {
		light: {
			palette: {
				primary: {
					main: "#1976d2",
					light: "#42a5f5",
				},
				warning: {
					main: "#ff9800",
				},
				success: {
					main: "#4caf50",
				},
				background: {
					default: "#f4f6f8",
					paper: "#f5f5f5",
				},
			},
		},
		dark: {
			palette: {
				primary: {
					main: "#90caf9",
					light: "#42a5f5",
				},
				warning: {
					main: "#ff9800",
				},
				success: {
					main: "#4caf50",
				},
				background: {
					default: "#121212",
					paper: "#1e1e1e",
				},
			},
		},
	},
	spacing: 8,
	shape: {
		borderRadius: 8,
	},
});

const AppContent: React.FC = () => {
	const { mode, setMode } = useColorScheme();
	const [initialized, setInitialized] = useState(false);

	// Initialize mode only once from localStorage or system preference
	useEffect(() => {
		if (!initialized) {
			const savedMode = localStorage.getItem("theme");
			if (savedMode && (savedMode === "light" || savedMode === "dark")) {
				setMode(savedMode);
			} else if (window.matchMedia("(prefers-color-scheme: dark)").matches) {
				setMode("dark");
			} else {
				setMode("light");
			}
			setInitialized(true);
		}
	}, [setMode, initialized]);

	// Sync changes back to localStorage (only when user changes mode)
	useEffect(() => {
		if (mode && initialized) {
			localStorage.setItem("theme", mode);
			document.documentElement.setAttribute("data-mui-color-scheme", mode);
			document.documentElement.setAttribute("data-bs-theme", mode); // Keep for Bootstrap compatibility
		}
	}, [mode, initialized]);

	return <App />;
};

export const AppWithTheme: React.FC = () => {
	return (
		<ThemeProvider theme={theme}>
			<CssBaseline enableColorScheme />
			<AppContent />
		</ThemeProvider>
	);
};
