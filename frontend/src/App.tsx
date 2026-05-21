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
import GroupsPage from "./pages/GroupsPage";
import MainPage from "./pages/landingPage";
import TemplateSelectionPage from "./pages/templateSelection";

function FilesystemRedirect() {
	const { ownerId, roomName } = useParams<{ ownerId: string; roomName: string }>();
	if (!ownerId || !roomName) return <Navigate to="/" replace />;
	return <Navigate to={`/rooms/${ownerId}/${roomName}?panel=filesystem`} replace />;
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
		path: "/groups",
		element: <GroupsPage />,
	},
	{
		path: "/rooms/:ownerId/:roomName/files",
		element: <FilesystemRedirect />,
	},
	{
		path: "/rooms/:ownerId/:roomName",
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
