import CssBaseline from "@mui/material/CssBaseline";
import { ThemeProvider, createTheme } from "@mui/material/styles";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { RouterProvider, createBrowserRouter } from "react-router-dom";
import FilesystemBrowserPage from "./pages/filesystemBrowser";
import MainPage from "./pages/landingPage";
import RoomListPage from "./pages/roomList";
import TemplateSelectionPage from "./pages/templateSelection";

const queryClient = new QueryClient({
	defaultOptions: {
		queries: {
			staleTime: 30000, // 30 seconds - matches existing hook usage
			gcTime: 5 * 60 * 1000, // 5 minutes
		},
	},
});

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
		path: "/rooms/:roomId/files",
		element: <FilesystemBrowserPage />,
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
