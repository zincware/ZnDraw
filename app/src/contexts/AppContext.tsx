import type React from "react";
import {
	type ReactNode,
	createContext,
	useContext,
	useEffect,
	useState,
} from "react";
import * as THREE from "three";
import type { Frame } from "../components/particles";
import type { IndicesState } from "../components/utils";
import type { HSLColor } from "../components/utils";
import { useTokenFromUrl } from "../hooks/useTokenFromUrl";
import { useVectorManager } from "../hooks/useVectorManager";
import type { RoomConfig } from "../types/room-config";
import { DEFAULT_ROOM_CONFIG } from "../types/room-config";

interface CameraAndControls {
	camera: THREE.Vector3;
	target: THREE.Vector3;
}

interface AppState {
	// Connection state
	connected: boolean;
	setConnected: (connected: boolean) => void;

	// Room state
	token: string;
	setToken: (token: string) => void;
	roomName: string;
	setRoomName: (roomName: string) => void;
	roomLock: boolean;
	setRoomLock: (roomLock: boolean) => void;

	// Frame state
	currentFrame: Frame;
	setCurrentFrame: (frame: Frame) => void;
	length: number;
	setLength: (length: number) => void;
	step: number;
	setStep: (step: number) => void;

	// Player state
	playing: boolean;
	setPlaying: (playing: boolean) => void;
	isFrameRendering: boolean;
	setIsFrameRendering: (rendering: boolean) => void;
	frameRate: number;
	setFrameRate: (rate: number) => void;

	// Selection state
	selectedFrames: IndicesState;
	setSelectedFrames: (frames: IndicesState) => void;
	selectedIds: Set<number>;
	setSelectedIds: (ids: Set<number>) => void;
	hoveredId: number;
	setHoveredId: (id: number | null) => void;

	// Drawing state
	isDrawing: boolean;
	setIsDrawing: (drawing: boolean) => void;
	points: THREE.Vector3[];
	setPoints: (points: THREE.Vector3[]) => void;
	selectedPoint: THREE.Vector3 | null;
	setSelectedPoint: (point: THREE.Vector3 | null) => void;

	// Camera state
	cameraAndControls: CameraAndControls;
	setCameraAndControls: (camera: CameraAndControls) => void;

	// Edit mode
	editMode: string;
	setEditMode: (mode: string) => void;

	// UI state
	colorMode: string;
	handleColorMode: () => void;
	showParticleInfo: boolean;
	setShowParticleInfo: (show: boolean) => void;
	tutorialURL: string;
	setTutorialURL: (url: string) => void;
	showSiMGen: boolean;
	setShowSiMGen: (show: boolean) => void;

	// Additional UI state
	modifierQueue: unknown[];
	isAuthenticated: boolean;
	setIsAuthenticated: (authenticated: boolean) => void;
	addPlotsWindow: number;
	setAddPlotsWindow: (show: number) => void;

	// Data collections
	bookmarks: Record<number, string>;
	setBookmarks: (bookmarks: Record<number, string>) => void;
	geometries: unknown[];
	setGeometries: (geometries: unknown[]) => void;

	// Configuration
	roomConfig: RoomConfig;
	setRoomConfig: (config: RoomConfig) => void;

	// Plotting
	updatedPlotsList: string[];
	setUpdatedPlotsList: (list: string[]) => void;

	// Messages
	messages: unknown[];
	setMessages: (messages: unknown[]) => void;

	// Vector data
	vectorProperties: Record<string, unknown>;
	setVectorProperties: (properties: Record<string, unknown>) => void;
	perParticleVectors: {
		start: THREE.Vector3;
		end: THREE.Vector3;
		vectorType: string;
	}[];
	vectorFieldData: [number, number, number][][];
	// Computed colormap for vectors
	vectorColormap: HSLColor[];
}

const AppContext = createContext<AppState | undefined>(undefined);

export const useAppContext = () => {
	const context = useContext(AppContext);
	if (context === undefined) {
		throw new Error("useAppContext must be used within an AppProvider");
	}
	return context;
};

interface AppProviderProps {
	children: ReactNode;
	colorMode: string;
	handleColorMode: () => void;
}

export const AppProvider: React.FC<AppProviderProps> = ({
	children,
	colorMode,
	handleColorMode,
}) => {
	// Extract token from URL
	const urlToken = useTokenFromUrl();

	// Connection state
	const [connected, setConnected] = useState<boolean>(false);

	// Room state - use urlToken as initial value and fallback to socket updates
	const [token, setToken] = useState<string>(urlToken || "");
	const [roomName, setRoomName] = useState<string>("");
	const [roomLock, setRoomLock] = useState<boolean>(false);

	// Frame state
	const [currentFrame, setCurrentFrame] = useState<Frame>({
		arrays: { colors: [], radii: [] },
		calc: null,
		cell: [],
		connectivity: [],
		info: null,
		numbers: [],
		pbc: [],
		positions: [],
		vectors: [],
		constraints: [],
	});
	const [length, setLength] = useState<number>(0);
	const [step, setStep] = useState<number>(0);

	// Player state
	const [playing, setPlaying] = useState<boolean>(false);
	const [isFrameRendering, setIsFrameRendering] = useState<boolean>(false);
	const [frameRate, setFrameRate] = useState<number>(1); // 1 = every frame, 2 = every second frame, etc.

	// Selection state
	const [selectedFrames, setSelectedFrames] = useState<IndicesState>({
		active: true,
		indices: new Set<number>(),
	});
	const [selectedIds, setSelectedIds] = useState<Set<number>>(new Set());
	const [hoveredId, setHoveredIdState] = useState<number>(-1);

	const setHoveredId = (id: number | null) => {
		setHoveredIdState(id ?? -1);
	};

	// Drawing state
	const [isDrawing, setIsDrawing] = useState<boolean>(false);
	const [points, setPoints] = useState<THREE.Vector3[]>([]);
	const [selectedPoint, setSelectedPoint] = useState<THREE.Vector3 | null>(
		null,
	);

	// Camera state
	const [cameraAndControls, setCameraAndControls] = useState<CameraAndControls>(
		{
			camera: new THREE.Vector3(5, 5, 5),
			target: new THREE.Vector3(0, 0, 0),
		},
	);

	// Edit mode
	const [editMode, setEditMode] = useState<string>("none");

	// UI state
	const [showParticleInfo, setShowParticleInfo] = useState<boolean>(false);
	const [tutorialURL, setTutorialURL] = useState<string>("");
	const [showSiMGen, setShowSiMGen] = useState<boolean>(false);

	// Additional UI state
	const [modifierQueue] = useState<unknown[]>([]);
	const [isAuthenticated, setIsAuthenticated] = useState<boolean>(false);
	const [addPlotsWindow, setAddPlotsWindow] = useState<number>(0);

	// Data collections
	const [bookmarks, setBookmarks] = useState<Record<number, string>>({});
	const [geometries, setGeometries] = useState<unknown[]>([]);

	// Configuration - using defaults from pydantic models
	const [roomConfig, setRoomConfig] = useState<RoomConfig>(DEFAULT_ROOM_CONFIG);

	// Enhanced config setter for debugging
	const setRoomConfigWithLogging = (config: RoomConfig) => {
		setRoomConfig(config);
	};

	// Plotting
	const [updatedPlotsList, setUpdatedPlotsList] = useState<string[]>([]);

	// Messages
	const [messages, setMessages] = useState<unknown[]>([]);

	// Vector data using custom hook
	const {
		vectorProperties,
		setVectorProperties,
		perParticleVectors,
		vectorFieldData,
		vectorColormap,
	} = useVectorManager({
		token,
		step,
		currentFrame,
		vectorConfig: roomConfig.VectorDisplay,
	});

	// Update token when urlToken changes from URL (prioritize URL token over socket token)
	useEffect(() => {
		if (urlToken && urlToken !== token) {
			setToken(urlToken);
			console.log("Token updated from URL:", urlToken);
		}
	}, [urlToken, token]);

	const value: AppState = {
		// Connection state
		connected,
		setConnected,

		// Room state
		token,
		setToken,
		roomName,
		setRoomName,
		roomLock,
		setRoomLock,

		// Frame state
		currentFrame,
		setCurrentFrame,
		length,
		setLength,
		step,
		setStep,

		// Player state
		playing,
		setPlaying,
		isFrameRendering,
		setIsFrameRendering,
		frameRate,
		setFrameRate,

		// Selection state
		selectedFrames,
		setSelectedFrames,
		selectedIds,
		setSelectedIds,
		hoveredId,
		setHoveredId,

		// Drawing state
		isDrawing,
		setIsDrawing,
		points,
		setPoints,
		selectedPoint,
		setSelectedPoint,

		// Camera state
		cameraAndControls,
		setCameraAndControls,

		// Edit mode
		editMode,
		setEditMode,

		// UI state
		colorMode,
		handleColorMode,
		showParticleInfo,
		setShowParticleInfo,
		tutorialURL,
		setTutorialURL,
		showSiMGen,
		setShowSiMGen,

		// Additional UI state
		modifierQueue,
		isAuthenticated,
		setIsAuthenticated,
		addPlotsWindow,
		setAddPlotsWindow,

		// Data collections
		bookmarks,
		setBookmarks,
		geometries,
		setGeometries,

		// Configuration
		roomConfig,
		setRoomConfig: setRoomConfigWithLogging,

		// Plotting
		updatedPlotsList,
		setUpdatedPlotsList,

		// Messages
		messages,
		setMessages,

		// Vector data
		vectorProperties,
		setVectorProperties,
		perParticleVectors,
		vectorFieldData,
		vectorColormap,
	};

	return <AppContext.Provider value={value}>{children}</AppContext.Provider>;
};
