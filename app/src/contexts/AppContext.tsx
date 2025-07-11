import React, { createContext, useContext, useState, ReactNode } from "react";
import * as THREE from "three";
import type { Frame } from "../components/particles";
import type { IndicesState } from "../components/utils";
import type { HSLColor } from "../components/utils";
import { useVectorManager } from "../hooks/useVectorManager";

// Define interfaces for better type safety
interface RoomConfig {
	Arrows: {
		[key: string]: any;
	};
	Particle: {
		[key: string]: any;
	};
	Visualization: {
		floor: boolean;
		[key: string]: any;
	};
	Camera: {
		camera_far: number;
		[key: string]: any;
	};
	PathTracer: {
		enabled: boolean;
		environment: string;
		[key: string]: any;
	};
	VectorDisplay: {
		vectorfield: boolean;
		vectors: string[];
		vector_scale: number;
		[key: string]: any;
	};
}

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

	// Selection state
	selectedFrames: IndicesState;
	setSelectedFrames: (frames: IndicesState) => void;
	selectedIds: Set<number>;
	setSelectedIds: (ids: Set<number>) => void;
	hoveredId: number;
	setHoveredId: (id: number) => void;

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
	modifierQueue: any[];
	isAuthenticated: boolean;
	setIsAuthenticated: (authenticated: boolean) => void;
	setAddPlotsWindow: (show: boolean) => void;

	// Data collections
	bookmarks: Record<number, string>;
	setBookmarks: (bookmarks: Record<number, string>) => void;
	geometries: any[];
	setGeometries: (geometries: any[]) => void;

	// Configuration
	roomConfig: RoomConfig;
	setRoomConfig: (config: RoomConfig) => void;

	// Plotting
	updatedPlotsList: string[];
	setUpdatedPlotsList: (list: string[]) => void;

	// Messages
	messages: any[];
	setMessages: (messages: any[]) => void;

	// Vector data
	vectorProperties: { [key: string]: any };
	setVectorProperties: (properties: { [key: string]: any }) => void;
	perParticleVectors: { start: THREE.Vector3; end: THREE.Vector3 }[];
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
	// Connection state
	const [connected, setConnected] = useState<boolean>(false);

	// Room state
	const [token, setToken] = useState<string>("");
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

	// Selection state
	const [selectedFrames, setSelectedFrames] = useState<IndicesState>({
		active: true,
		indices: new Set<number>(),
	});
	const [selectedIds, setSelectedIds] = useState<Set<number>>(new Set());
	const [hoveredId, setHoveredId] = useState<number>(-1);

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
	const [modifierQueue] = useState<any[]>([]);
	const [isAuthenticated, setIsAuthenticated] = useState<boolean>(false);
	const [addPlotsWindow, setAddPlotsWindow] = useState<boolean>(false);

	// Data collections
	const [bookmarks, setBookmarks] = useState<Record<number, string>>({});
	const [geometries, setGeometries] = useState<any[]>([]);

	// Configuration
	const [roomConfig, setRoomConfig] = useState<RoomConfig>({
		Arrows: {
			colormap: "viridis",
			colorrange: [0, 1],
			normalize: false,
			opacity: 1.0,
			rescale: 1.0,
			scale_vector_thickness: false,
		},
		Particle: {
			bond_size: 1.0,
			hover_opacity: 0.6,
			material: "MeshStandardMaterial",
			particle_size: 1.0,
			selection_color: "#ff5722",
			selection_opacity: 0.7,
		},
		Visualization: {
			camera_target: null,
			floor: false,
			frame_update: true,
			fps: 24,
			loop: true,
			material: "MeshStandardMaterial",
		},
		Camera: {
			camera_far: 100,
			camera_fov: 75,
			camera_near: 0.1,
			controls: "OrbitControls",
			position: [10, 10, 10],
			synchronize_camera: false,
			up: [0, 1, 0],
		},
		PathTracer: {
			enabled: false,
			environment: "city",
		},
		VectorDisplay: {
			colormap: "viridis",
			colorrange: [0, 1],
			normalize: false,
			opacity: 1.0,
			rescale: 1.0,
			scale_vector_thickness: false,
			vectorfield: false,
			vectors: [],
			vector_scale: 1.0,
		},
	});

	// Plotting
	const [updatedPlotsList, setUpdatedPlotsList] = useState<string[]>([]);

	// Messages
	const [messages, setMessages] = useState<any[]>([]);

	// Vector data using custom hook
	const { vectorProperties, setVectorProperties, perParticleVectors, vectorFieldData, vectorColormap } = useVectorManager({
		token,
		step,
		currentFrame,
		vectorConfig: roomConfig.VectorDisplay,
	});

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
		setAddPlotsWindow,

		// Data collections
		bookmarks,
		setBookmarks,
		geometries,
		setGeometries,

		// Configuration
		roomConfig,
		setRoomConfig,

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
