import { useAppStore } from '../store';
// import { usePlaybackControls } from '../hooks/usePlaybackControls'; // This hook would manage the presenter token

export function TimelineScrubber() {
  const { currentFrame, frameCount } = useAppStore(state => ({
    currentFrame: state.currentFrame,
    frameCount: state.frameCount,
  }));
  
  // This custom hook hides all the socket.emit complexity
  // const { changeFrame, startScrubbing, stopScrubbing } = usePlaybackControls();

  // const handleScrub = (event) => {
  //   const newFrame = parseInt(event.target.value, 10);
  //   // The hook checks if we are the presenter before sending the update
  //   changeFrame(newFrame);
  // };

  return (
    <div className="scrubber-container">
      <span>{currentFrame}</span>
      <input
        type="range"
        min="0"
        max={frameCount > 0 ? frameCount - 1 : 0}
        value={currentFrame}
        // onMouseDown={startScrubbing} // Emits 'request_presenter_token'
        // onMouseUp={stopScrubbing}   // Emits 'release_presenter_token'
        // onChange={handleScrub}      // Emits 'set_room_frame'
      />
      <span>{frameCount}</span>
    </div>
  );
}