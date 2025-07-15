import { CssBaseline } from "@mui/material";
import { ThemeProvider, createTheme } from "@mui/material/styles";
import React from "react";
import ReactDOM from "react-dom/client";
import App from "./App.tsx";
import "./index.css";

const theme = createTheme({
  palette: {
    primary: {
      main: '#1976d2', // Example primary color
      light: '#42a5f5',
    },
    warning: {
      main: '#ff9800', // Example warning color for bookmarks
    },
    success: {
        main: '#4caf50', // Example success color for Wifi icon
    },
    background: {
      default: '#f4f6f8', // Example background color
      paper: '#ffffff',
    },
  },
  // You can define spacing, typography, shadows, etc. here
  spacing: 8, // Default spacing unit (8px)
  shape: {
    borderRadius: 8, // Default border radius
  },
  shadows: [
    'none', // shadow[0]
    '0px 2px 1px -1px rgba(0,0,0,0.2),0px 1px 1px 0px rgba(0,0,0,0.14),0px 1px 3px 0px rgba(0,0,0,0.12)', // shadow[1]
    // ... and so on for other shadows
  ],
});


// React strict mode renders the app twice to detect side effects
// this will fail for our useRef based socket detection
// and messages will be send through the socket, that should not be send
ReactDOM.createRoot(document.getElementById("root")!).render(
	<React.StrictMode>
		<ThemeProvider theme={theme}>
			<CssBaseline />
			<App />
		</ThemeProvider>
	</React.StrictMode>,
);
