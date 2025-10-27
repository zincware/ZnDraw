import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import Box from "@mui/material/Box";
import Container from "@mui/material/Container";
import Typography from "@mui/material/Typography";
import CircularProgress from "@mui/material/CircularProgress";
import Alert from "@mui/material/Alert";
import { listRooms, getDefaultRoom } from "../myapi/client";

/**
 * StartupPage implements the new startup logic:
 * - No rooms: Navigate to empty template
 * - One room: Navigate to that room
 * - Multiple rooms with default: Navigate to default room
 * - Multiple rooms without default: Show room list
 */
export default function TemplateSelectionPage() {
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const navigate = useNavigate();

  useEffect(() => {
    const determineStartupNavigation = async () => {
      try {
        const rooms = await listRooms();

        if (rooms.length === 0) {
          // No rooms - create empty template
          const roomUuid = crypto.randomUUID();
          navigate(`/rooms/${roomUuid}?template=empty`);
        } else if (rooms.length === 1) {
          // One room - navigate to it
          navigate(`/rooms/${rooms[0].id}`);
        } else {
          // Multiple rooms - check for default
          const { roomId: defaultRoomId } = await getDefaultRoom();
          
          if (defaultRoomId) {
            // Navigate to default room
            navigate(`/rooms/${defaultRoomId}`);
          } else {
            // No default - show room list
            navigate("/rooms");
          }
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
            Loading rooms...
          </Typography>
        </Box>
      </Container>
    );
  }

  if (error) {
    return (
      <Container maxWidth="md">
        <Box sx={{ mt: 4 }}>
          <Alert severity="error">
            Failed to load rooms: {error}
          </Alert>
        </Box>
      </Container>
    );
  }

  // This component only handles navigation logic, so if we reach here
  // something went wrong with the navigation
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
