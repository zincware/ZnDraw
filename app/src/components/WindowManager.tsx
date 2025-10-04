import { useWindowManagerStore } from "../stores/windowManagerStore";
import FigureWindow from "./FigureWindow";
import { useShallow } from "zustand/react/shallow";

function WindowManager() {
  const windowIds = useWindowManagerStore(
    useShallow((state) => Object.keys(state.openWindows)),
  );

  return (
    <>
      {windowIds.map((windowId) => (
        <FigureWindow key={windowId} windowId={windowId} />
      ))}
    </>
  );
}

export default WindowManager;
