import {
  createBrowserRouter,
  RouterProvider,
  useSearchParams,
} from "react-router-dom";
import MainPage from "./pages/landingPage";
import TemplateSelectionPage from "./pages/templateSelection";
import RoomListPage from "./pages/roomList";
import RoomWaitingPage from "./pages/roomWaiting";
import FileBrowserPage from "./pages/fileBrowser";
import RemoteFileBrowserPage from "./pages/remoteFileBrowser";
import { ExtensionsOverview } from "./pages/extensionsOverview";
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

// Component for /rooms/:roomId route
// Shows waiting page if waitForCreation=true, otherwise shows main page
const RoomPage = () => {
  const [searchParams] = useSearchParams();
  const waitForCreation = searchParams.get("waitForCreation") === "true";

  if (waitForCreation) {
    return <RoomWaitingPage />;
  }

  return <MainPage />;
};

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
    path: "/extensions",
    element: <ExtensionsOverview mode="global" />,
  },
  {
    path: "/rooms/:roomId/extensions",
    element: <ExtensionsOverview mode="room" />,
  },
  {
    path: "/rooms/:roomId",
    element: <RoomPage />,
  },
  {
    path: "/room/:roomId",
    element: <RoomPage />,
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
