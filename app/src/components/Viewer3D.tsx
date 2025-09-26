import { useEffect } from 'react';
import { useAppStore } from '../store';
import { useTrajectoryData } from '../hooks/useTrajectoryData';

export function Viewer3D() {
  const currentFrame = useAppStore(state => state.currentFrame);
  const { loadFrame, frameDataCache } = useTrajectoryData();

  // Effect to load data when the current frame changes
  useEffect(() => {
    // Define what data this component needs for the current frame
    const requiredKeys = ['positions', 'species'];
    loadFrame(currentFrame, requiredKeys);
  }, [currentFrame, loadFrame]);

  const data = frameDataCache.get(currentFrame);

  if (!data || !data.positions) {
    return <div>Loading frame {currentFrame}...</div>;
  }

  // Your actual 3D rendering logic goes here
  return (
    <div className="viewer">
      {/* e.g., <ThreeCanvas positions={data.positions} species={data.species} /> */}
      <p>Displaying positions for frame {currentFrame}.</p>
    </div>
  );
}