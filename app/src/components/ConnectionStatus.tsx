import { useAppStore } from "../store";

export function ConnectionStatus() {
  // Use individual selectors to prevent unnecessary re-renders
  const isConnected = useAppStore((state) => state.isConnected);
  const currentFrame = useAppStore((state) => state.currentFrame);
  const frameCount = useAppStore((state) => state.frameCount);
  const roomId = useAppStore((state) => state.roomId);

  return (
    <div className="status">
      <p>
        Status: {isConnected ? "Connected" : "Disconnected"}
        {isConnected && <span> - Total Frames: {frameCount}</span>}
        {isConnected && <span> - Current Frame: {currentFrame}</span>}
        {roomId && <span> - Room: {roomId}</span>}
      </p>
    </div>
  );
}
