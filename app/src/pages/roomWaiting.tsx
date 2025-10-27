import { useEffect, useState } from "react";
import { useParams, useNavigate, Navigate } from "react-router-dom";
import {
  Box,
  CircularProgress,
  Typography,
  Alert,
  Button,
} from "@mui/material";
import { getRoomInfo } from "../myapi/client";

const POLL_INTERVAL_MS = 500; // Poll every 0.5 seconds
const TIMEOUT_MS = 60000; // Timeout after 60 seconds

export default function RoomWaitingPage() {
  const { roomId } = useParams<{ roomId: string }>();
  const navigate = useNavigate();
  const [status, setStatus] = useState<"waiting" | "ready" | "timeout" | "error">("waiting");
  const [errorMessage, setErrorMessage] = useState<string>("");
  const [elapsedTime, setElapsedTime] = useState(0);

  useEffect(() => {
    if (!roomId) {
      setStatus("error");
      setErrorMessage("No room ID provided");
      return;
    }

    let pollTimer: NodeJS.Timeout;
    let elapsedTimer: NodeJS.Timeout;
    let isActive = true;
    const startTime = Date.now();

    const checkRoomStatus = async () => {
      try {
        const roomInfo = await getRoomInfo(roomId);
        
        // Room exists and has frames - ready to join!
        if (roomInfo.frameCount > 0) {
          if (isActive) {
            setStatus("ready");
            // Small delay to show success state before redirect
            setTimeout(() => {
              navigate(`/rooms/${roomId}`, { replace: true });
            }, 500);
          }
          return;
        }

        // Room exists but no frames yet - keep waiting
        if (isActive) {
          const elapsed = Date.now() - startTime;
          if (elapsed >= TIMEOUT_MS) {
            setStatus("timeout");
          } else {
            // Continue polling
            pollTimer = setTimeout(checkRoomStatus, POLL_INTERVAL_MS);
          }
        }
      } catch (error: any) {
        // 404 means room doesn't exist yet - keep waiting
        if (error?.response?.status === 404) {
          if (isActive) {
            const elapsed = Date.now() - startTime;
            if (elapsed >= TIMEOUT_MS) {
              setStatus("timeout");
            } else {
              // Continue polling
              pollTimer = setTimeout(checkRoomStatus, POLL_INTERVAL_MS);
            }
          }
        } else {
          // Other errors
          if (isActive) {
            setStatus("error");
            setErrorMessage(
              error?.response?.data?.message || 
              error?.message || 
              "Failed to check room status"
            );
          }
        }
      }
    };

    // Start polling
    checkRoomStatus();

    // Update elapsed time display every second
    elapsedTimer = setInterval(() => {
      setElapsedTime(Math.floor((Date.now() - startTime) / 1000));
    }, 1000);

    // Cleanup
    return () => {
      isActive = false;
      clearTimeout(pollTimer);
      clearInterval(elapsedTimer);
    };
  }, [roomId, navigate]);

  if (!roomId) {
    return <Navigate to="/rooms" replace />;
  }

  return (
    <Box
      sx={{
        display: "flex",
        flexDirection: "column",
        alignItems: "center",
        justifyContent: "center",
        height: "100vh",
        width: "100vw",
        padding: 4,
        gap: 3,
      }}
    >
      {status === "waiting" && (
        <>
          <CircularProgress size={60} />
          <Typography variant="h5" align="center">
            Loading file into room...
          </Typography>
          <Typography variant="body1" color="text.secondary" align="center">
            Room: <strong>{roomId}</strong>
          </Typography>
          <Typography variant="body2" color="text.secondary">
            Elapsed: {elapsedTime}s
          </Typography>
          <Typography variant="caption" color="text.secondary" align="center">
            This usually takes a few seconds. Large files may take longer.
          </Typography>
        </>
      )}

      {status === "ready" && (
        <>
          <CircularProgress size={60} color="success" />
          <Typography variant="h5" color="success.main" align="center">
            Room ready! Redirecting...
          </Typography>
        </>
      )}

      {status === "timeout" && (
        <>
          <Alert severity="warning" sx={{ maxWidth: 600 }}>
            <Typography variant="h6" gutterBottom>
              Timeout
            </Typography>
            <Typography variant="body2" paragraph>
              The room <strong>{roomId}</strong> is taking longer than expected to load.
              This could mean:
            </Typography>
            <ul style={{ marginTop: 8, marginBottom: 8 }}>
              <li>The file is very large and still loading</li>
              <li>There was an issue with the file upload</li>
              <li>The Celery worker may not be running</li>
            </ul>
            <Typography variant="body2">
              You can try waiting longer or check the server logs.
            </Typography>
          </Alert>
          <Box sx={{ display: "flex", gap: 2 }}>
            <Button
              variant="contained"
              onClick={() => {
                navigate(`/rooms/${roomId}`, { replace: true });
              }}
            >
              Try to Join Anyway
            </Button>
            <Button
              variant="outlined"
              onClick={() => navigate("/rooms", { replace: true })}
            >
              Back to Room List
            </Button>
          </Box>
        </>
      )}

      {status === "error" && (
        <>
          <Alert severity="error" sx={{ maxWidth: 600 }}>
            <Typography variant="h6" gutterBottom>
              Error
            </Typography>
            <Typography variant="body2">
              {errorMessage || "An unexpected error occurred"}
            </Typography>
          </Alert>
          <Button
            variant="contained"
            onClick={() => navigate("/rooms", { replace: true })}
          >
            Back to Room List
          </Button>
        </>
      )}
    </Box>
  );
}
