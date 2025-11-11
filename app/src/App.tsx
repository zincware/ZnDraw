import {
  createBrowserRouter,
  RouterProvider,
} from "react-router-dom";
import MainPage from "./pages/landingPage";
import TemplateSelectionPage from "./pages/templateSelection";
import RoomListPage from "./pages/roomList";
import FileBrowserPage from "./pages/fileBrowser";
import RemoteFileBrowserPage from "./pages/remoteFileBrowser";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ThemeProvider, createTheme } from "@mui/material/styles";
import CssBaseline from "@mui/material/CssBaseline";

const queryClient = new QueryClient();

const theme = createTheme({
  colorSchemes: {
    dark: true,
  },
  defaultColorScheme: "light",
});

const router = createBrowserRouter([
  {
    path: "/",
    element: <TemplateSelectionPage />,
  },
  {
    path: "/rooms",
    element: <RoomListPage />,
  },
  {
    path: "/file-browser",
    element: <FileBrowserPage />,
  },
  {
    path: "/rooms/:roomId/remote-files",
    element: <RemoteFileBrowserPage />,
  },
  {
    path: "/rooms/:roomId",
    element: <MainPage />,
  },
  {
    path: "/room/:roomId",
    element: <MainPage />,
  },
]);

export function App() {
  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <QueryClientProvider client={queryClient}>
        <RouterProvider router={router} />
      </QueryClientProvider>
    </ThemeProvider>
  );
}
