import { useSocketManager } from '../hooks/useSocketManager';
import { ConnectionStatus } from '../components/ConnectionStatus';
import { useTrajectoryData } from '../hooks/useTrajectoryData';
import { useParams } from 'react-router-dom'; // 1. Import useParams


export function RoomPage() {
  // Initialize the socket manager once. It works in the background.
  const { roomId } = useParams<{ roomId: string }>();
  useSocketManager(roomId as string);
  const { loadFrame } = useTrajectoryData();

  return (
    <div className="app-layout">
      <header>
        <h1>Collaborative Viewer</h1>
        <ConnectionStatus />
        <button onClick={() => loadFrame(0, ["positions", "numbers"])}>Load Frame 0</button>
      </header>
    </div>
  );
}

export default RoomPage;