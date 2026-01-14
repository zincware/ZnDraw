import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import Box from "@mui/material/Box";
import Container from "@mui/material/Container";
import Typography from "@mui/material/Typography";
import CircularProgress from "@mui/material/CircularProgress";
import Alert from "@mui/material/Alert";
import {
	listRooms,
	getDefaultRoom,
	getRoom,
	duplicateRoom,
} from "../myapi/client";
import { getLastVisitedRoom } from "../utils/roomTracking";
import { ensureAuthenticated } from "../utils/auth";

/**
 * Check if a room exists by attempting to fetch it.
 *
 * Parameters
 * ----------
 * roomId : string
 *     The room ID to check.
 *
 * Returns
 * -------
 * boolean
 *     True if room exists, false otherwise.
 */
async function roomExists(roomId: string): Promise<boolean> {
	try {
		await getRoom(roomId);
		return true;
	} catch {
		return false;
	}
}

/**
 * Determine which room to use as template.
 *
 * Auto-template: If only 1 room exists with data, use it.
 * Otherwise: Use explicit template from getDefaultRoom().
 *
 * Returns
 * -------
 * string | null
 *     Template room ID or null if no template.
 */
async function determineTemplate(): Promise<string | null> {
	const rooms = await listRooms();

	// Auto-template: single room with data
	if (rooms.length === 1 && rooms[0].frameCount > 0) {
		console.log("[Template] Auto-template from single room:", rooms[0].id);
		return rooms[0].id;
	}

	// Explicit template
	const { roomId: explicitTemplate } = await getDefaultRoom();
	if (explicitTemplate) {
		const exists = await roomExists(explicitTemplate);
		if (exists) {
			console.log("[Template] Using explicit template:", explicitTemplate);
			return explicitTemplate;
		}
	}

	console.log("[Template] No template found");
	return null;
}

/**
 * TemplateSelectionPage implements simple startup navigation:
 *
 * 1. Check localStorage for last visited room
 *    - If exists and valid -> navigate to it
 *    - Otherwise -> create new room
 *
 * 2. Create new room:
 *    - Determine template (auto or explicit)
 *    - Duplicate from template or create empty
 *    - Navigate to new room
 */
export default function TemplateSelectionPage() {
	const [loading, setLoading] = useState(true);
	const [error, setError] = useState<string | null>(null);
	const navigate = useNavigate();

	useEffect(() => {
		const determineStartupNavigation = async () => {
			try {
				// Ensure user is authenticated before making any API calls
				await ensureAuthenticated();

				// Step 1: Check localStorage for last visited room
				const lastRoomId = getLastVisitedRoom();
				console.log("[Startup] Last room from localStorage:", lastRoomId);

				if (lastRoomId) {
					const exists = await roomExists(lastRoomId);
					console.log("[Startup] Last room exists:", exists);
					if (exists) {
						console.log("[Startup] Navigating to last room:", lastRoomId);
						navigate(`/rooms/${lastRoomId}`);
						return;
					}
					// Room was deleted, continue to create new
					console.log("[Startup] Last room deleted, creating new room");
				}

				// Step 2: Create new room (from template or empty)
				const templateRoomId = await determineTemplate();
				const newRoomId = crypto.randomUUID();

				if (templateRoomId) {
					console.log(
						"[Startup] Creating room from template:",
						templateRoomId,
						"->",
						newRoomId,
					);
					await duplicateRoom(templateRoomId, { newRoomId });
					navigate(`/rooms/${newRoomId}`);
				} else {
					console.log("[Startup] Creating default room:", newRoomId);
					navigate(`/rooms/${newRoomId}`);
				}
			} catch (err) {
				setError(err instanceof Error ? err.message : "Unknown error");
			} finally {
				setLoading(false);
			}
		};

		determineStartupNavigation();
	}, [navigate]);

	if (loading) {
		return (
			<Container maxWidth="md">
				<Box
					sx={{
						display: "flex",
						justifyContent: "center",
						alignItems: "center",
						minHeight: "100vh",
					}}
				>
					<CircularProgress />
					<Typography variant="body1" sx={{ ml: 2 }}>
						Loading...
					</Typography>
				</Box>
			</Container>
		);
	}

	if (error) {
		return (
			<Container maxWidth="md">
				<Box sx={{ mt: 4 }}>
					<Alert severity="error">Failed to load: {error}</Alert>
				</Box>
			</Container>
		);
	}

	// This component only handles navigation logic
	return (
		<Container maxWidth="md">
			<Box sx={{ mt: 4 }}>
				<Alert severity="warning">
					Navigation logic did not redirect properly. Please refresh the page.
				</Alert>
			</Box>
		</Container>
	);
}
