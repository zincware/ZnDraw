import { useWindowManagerStore } from '../stores/windowManagerStore';
import FigureWindow from './FigureWindow';

function WindowManager() {
  const { openWindows } = useWindowManagerStore();

  return (
    <>
      {Object.keys(openWindows).map((figureKey) => (
        <FigureWindow key={figureKey} figureKey={figureKey} />
      ))}
    </>
  );
}

export default WindowManager;