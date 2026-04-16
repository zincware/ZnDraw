import CssBaseline from "@mui/material/CssBaseline";
import { ThemeProvider, createTheme } from "@mui/material/styles";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import {
	Navigate,
	RouterProvider,
	createBrowserRouter,
	useParams,
} from "react-router-dom";
import { MuiCssVars } from "./MuiCssVars";
import CliLoginApprovePage from "./pages/cliLoginApprove";
import MainPage from "./pages/landingPage";
import RoomListPage from "./pages/roomList";
import TemplateSelectionPage from "./pages/templateSelection";

function FilesystemRedirect() {
	const { roomId } = useParams<{ roomId: string }>();
	if (!roomId) return <Navigate to="/" replace />;
	return <Navigate to={`/rooms/${roomId}?panel=filesystem`} replace />;
}

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
		path: "/auth/cli",
		element: <CliLoginApprovePage />,
	},
	{
		path: "/rooms",
		element: <RoomListPage />,
	},
	{
		path: "/rooms/:roomId/files",
		element: <FilesystemRedirect />,
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
			<MuiCssVars />
			<QueryClientProvider client={queryClient}>
				<RouterProvider router={router} />
			</QueryClientProvider>
		</ThemeProvider>
	);
}
