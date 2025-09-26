import { useAppStore } from '../store';

export function ConnectionStatus() {
  const { isConnected, currentFrame, frameData, frameCount, roomId } = useAppStore();

  return (
    <div className="status">
      <p>
        Status: {isConnected ? 'Connected' : 'Disconnected'}
        {isConnected && <span> - Total Frames: {frameCount}</span>}
        {isConnected && <span> - Current Frame: {currentFrame}</span>}
        {frameData && <span> - Frame Data: {JSON.stringify(frameData)}</span>}
        {roomId && <span> - Room: {roomId}</span>}
      </p>
    </div>
  );
}
