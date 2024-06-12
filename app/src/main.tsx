import React from "react";
import ReactDOM from "react-dom/client";
import App from "./App.tsx";
import "./index.css";

import "bootstrap/dist/css/bootstrap.min.css";

// React strict mode renders the app twice to detect side effects
// this will fail for our useRef based socket detection
// and messages will be send through the socket, that should not be send
ReactDOM.createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>,
);
