import { useEffect } from "react";
import {
	setupBookmarks,
	setupCamera,
	setupConfig,
	setupFigures,
	setupFrames,
	setupGeometries,
	setupMessages,
	setupPoints,
	setupSelection,
	setupStep,
} from "../components/api";
import { useAppContext } from "../contexts/AppContext";
import { socket } from "../socket";

export const useSocketManager = () => {
	const {
		connected,
		setConnected,
		setRoomLock,
		setToken,
		setTutorialURL,
		setShowSiMGen,
		isAuthenticated,
		setIsAuthenticated,
		// Frame and player state
		step,
		setStep,
		currentFrame,
		setCurrentFrame,
		length,
		setLength,
		playing,
		setPlaying,
		setIsFrameRendering,
		// Selection and interaction
		selectedIds,
		setSelectedIds,
		selectedFrames,
		// Data collections
		bookmarks,
		setBookmarks,
		points,
		setPoints,
		geometries,
		setGeometries,
		// Camera and config
		cameraAndControls,
		setCameraAndControls,
		roomConfig,
		setRoomConfig,
		// Messages and plotting
		messages,
		setMessages,
		setUpdatedPlotsList,
		// Token for API calls
		token,
	} = useAppContext();

	// Socket connection management
	useEffect(() => {
		function onConnect() {
			setConnected(true);
			console.log("connected");

			// Connect to room and get authentication status
			socket.emit(
				"webclient:connect",
				(data: { name: string; room: string; authenticated: boolean }) => {
					setIsAuthenticated(data.authenticated);
					setToken(data.room);
				},
			);

			// get lock state
			socket.emit("room:lock:get", (data: boolean) => {
				setRoomLock(data);
			});
		}

		function onDisconnect() {
			setConnected(false);
			setIsAuthenticated(false);
		}

		function onTutorialURL(data: string) {
			setTutorialURL(data);
		}

		function onShowSiMGen(data: boolean) {
			setShowSiMGen(data);
		}

		function onRoomLockSet(locked: boolean) {
			setRoomLock(locked);
		}

		socket.on("connect", onConnect);
		socket.on("disconnect", onDisconnect);
		socket.on("tutorial:url", onTutorialURL);
		socket.on("showSiMGen", onShowSiMGen);
		socket.on("room:lock:set", onRoomLockSet);

		return () => {
			socket.off("connect", onConnect);
			socket.off("disconnect", onDisconnect);
			socket.off("tutorial:url", onTutorialURL);
			socket.off("showSiMGen", onShowSiMGen);
			socket.off("room:lock:set", onRoomLockSet);
		};
	}, [
		setConnected,
		setIsAuthenticated,
		setRoomLock,
		setToken,
		setTutorialURL,
		setShowSiMGen,
	]);

	// Setup all the socket-based data management - these are custom hooks
	// They must be called at the top level, not inside useEffect
	setupStep(token, setStep, step);
	setupFrames(
		token,
		step,
		setCurrentFrame,
		currentFrame,
		setLength,
		setStep,
		roomConfig.Visualization.frame_update || true,
		setIsFrameRendering,
	);
	setupSelection(token, setSelectedIds, selectedIds);
	setupBookmarks(token, setBookmarks, bookmarks);
	setupPoints(token, setPoints, points);
	setupGeometries(token, setGeometries, geometries);
	setupCamera(
		token,
		cameraAndControls,
		setCameraAndControls,
		roomConfig.Camera.synchronize_camera || true,
	);
	setupConfig(token, setRoomConfig);
	setupFigures(token, setUpdatedPlotsList);
	setupMessages(token, setMessages, messages);

	return {
		connected,
		token,
	};
};
