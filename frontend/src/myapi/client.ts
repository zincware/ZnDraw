import axios from "axios";
import { useAppStore } from "../store";
import { acquireToken, getToken, logout } from "../utils/auth";
import { packBinary, unpackBinary } from "../utils/msgpack-numpy";
import type { ChatMessage, ChatMessagesResponse } from "../types/chat";
import type { TaskStatus } from "../types/jobs";

/**
 * Decodes frame data from the backend msgpack format.
 *
 * Backend sends: List[Dict[bytes, bytes]] where:
 * - Outer layer is msgpack list with bytes keys
 * - Each frame is dict with bytes keys and msgpack-numpy encoded bytes values
 * - Keys like b"arrays.positions", b"arrays.numbers", etc.
 *
 * Our custom unpackBinary handles:
 * - Bytes keys at all levels (using mapKeyConverter)
 * - Nested msgpack-numpy encoding (recursive decoding)
 * - Row-major array format (no transposition needed)
 *
 * @param encoded - Raw arraybuffer from HTTP response
 * @param key - Key to extract from frame (e.g., "arrays.positions")
 * @returns Decoded value for the specified key (TypedArray or other data)
 */
function decodeTypedData(encoded: any, key: string) {
	if (!encoded) return undefined;
	try {
		// Decode everything in one go with our custom msgpack-numpy decoder
		// This handles bytes keys, nested msgpack, and numpy arrays automatically
		const framesList = unpackBinary(
			encoded instanceof Uint8Array ? encoded : new Uint8Array(encoded),
		);

		if (!Array.isArray(framesList) || framesList.length === 0) {
			console.error("Expected non-empty array of frames");
			return undefined;
		}

		// Get first frame (assuming single frame query based on original code)
		const frameDict = framesList[0];

		if (typeof frameDict !== "object" || frameDict === null) {
			console.error("Expected frame to be a dict");
			return undefined;
		}

		// Get the value for the requested key
		// Keys are strings (converted from bytes by mapKeyConverter)
		const value = frameDict[key];

		if (value === undefined) {
			console.error(
				`Key "${key}" not found in frame. Available keys:`,
				Object.keys(frameDict),
			);
			return undefined;
		}

		// Value is already decoded (TypedArray for numpy arrays, or other types)
		return value;
	} catch (err) {
		console.error(`Failed to decode ${key}:`, err);
		return undefined;
	}
}

// Define the types based on your Python code
export interface FigureData {
	type: "plotly";
	[key: string]: any; // Allows for flexible figure data
}

export interface FigureResponse {
	key: string;
	figure: FigureData;
}

export interface FigureListResponse {
	items: string[];
}

export interface GlobalSettings {
	simgen: {
		enabled: boolean;
	};
}

const apiClient = axios.create({});

// Add interceptor to include JWT token, session ID, and lock token in all requests
apiClient.interceptors.request.use((config) => {
	const token = getToken();
	if (token) {
		config.headers["Authorization"] = `Bearer ${token}`;
	}

	const sessionId = useAppStore.getState().sessionId;
	if (sessionId) {
		config.headers["X-Session-ID"] = sessionId;
	}

	// Include Lock-Token for mutation requests when holding a lock
	const lockToken = useAppStore.getState().lockToken;
	if (lockToken) {
		config.headers["Lock-Token"] = lockToken;
	}

	return config;
});

// 401 response interceptor: acquire fresh token and retry once.
// acquireToken() coalesces concurrent calls internally.
const retriedRequests = new WeakSet<object>();

apiClient.interceptors.response.use(
	(response) => response,
	async (error) => {
		const originalRequest = error.config;

		// Only handle 401, skip if already retried or auth endpoint (avoid loops)
		if (
			error.response?.status !== 401 ||
			retriedRequests.has(originalRequest) ||
			originalRequest.url?.startsWith("/v1/auth/")
		) {
			return Promise.reject(error);
		}

		retriedRequests.add(originalRequest);

		// Clear the known-bad token so acquireToken() goes straight to guest login
		// instead of making a redundant validation call to /me.
		logout();
		const { token } = await acquireToken();
		useAppStore
			.getState()
			.showSnackbar("Session expired — reconnected as guest", "warning");
		originalRequest.headers["Authorization"] = `Bearer ${token}`;
		return apiClient(originalRequest);
	},
);

// --- API Functions ---

export const getServerVersion = async (): Promise<{ version: string }> => {
	const { data } = await apiClient.get("/v1/version");
	return data;
};

export const getGlobalSettings = async (): Promise<GlobalSettings> => {
	const { data } = await apiClient.get("/v1/config/global-settings");
	return data;
};

export const listFigures = async (
	roomId: string,
): Promise<FigureListResponse> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}/figures`);
	return data;
};

export const getFigure = async (
	roomId: string,
	key: string,
): Promise<FigureResponse> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}/figures/${key}`);
	return data;
};

export const createFigure = async (
	roomId: string,
	key: string,
	figure: FigureData,
): Promise<{ key: string; created: boolean }> => {
	const { data } = await apiClient.post(
		`/v1/rooms/${roomId}/figures/${encodeURIComponent(key)}`,
		{ figure },
	);
	return data;
};

export const deleteFigure = async (
	roomId: string,
	key: string,
): Promise<{ status: string }> => {
	const { data } = await apiClient.delete(`/v1/rooms/${roomId}/figures/${key}`);
	return data;
};

export interface PropertyMeta {
	dtype: string;
	shape?: number[];
	type: "array" | "scalar";
}

export interface FrameMetadata {
	frame_id: number;
	metadata: Record<string, PropertyMeta>;
	source_room: string;
}

export const getFrameMetadata = async (
	roomId: string,
	frameId = 0,
): Promise<FrameMetadata> => {
	const { data } = await apiClient.get(
		`/v1/rooms/${roomId}/frames/${frameId}/metadata`,
	);
	return data;
};

// ==================== Geometries API ====================

export interface GeometryData {
	type: string;
	data: Record<string, any>;
	selection: number[];
}

export interface GeometryResponse {
	key: string;
	geometry: GeometryData;
}

export interface GeometryListResponse {
	items: Record<string, GeometryData>;
	types?: {
		schemas: Record<string, any>;
		defaults: Record<string, any>;
	};
}

export const listGeometries = async (
	roomId: string,
): Promise<GeometryListResponse> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}/geometries`);
	return data;
};

export const getGeometry = async (
	roomId: string,
	key: string,
): Promise<GeometryResponse> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}/geometries/${key}`);
	return data;
};

export const createGeometry = async (
	roomId: string,
	key: string,
	geometryType: string,
	geometryData: Record<string, any>,
): Promise<{ status: string }> => {
	try {
		const { data } = await apiClient.put(
			`/v1/rooms/${roomId}/geometries/${key}`,
			{
				type: geometryType,
				data: geometryData,
			},
		);
		return data;
	} catch (error: any) {
		if (error.response?.status === 423) {
			useAppStore
				.getState()
				.showSnackbar(
					"Cannot modify geometry - room is locked by another user",
					"warning",
				);
		}
		throw error;
	}
};

export const deleteGeometry = async (
	roomId: string,
	key: string,
): Promise<{ status: string }> => {
	try {
		const { data } = await apiClient.delete(
			`/v1/rooms/${roomId}/geometries/${key}`,
		);
		return data;
	} catch (error: any) {
		if (error.response?.status === 423) {
			useAppStore
				.getState()
				.showSnackbar(
					"Cannot delete geometry - resource is locked by another user",
					"warning",
				);
		}
		throw error;
	}
};

// ==================== Default Camera API ====================

export interface DefaultCameraResponse {
	default_camera: string | null;
}

export const getDefaultCamera = async (
	roomId: string,
): Promise<DefaultCameraResponse> => {
	const { data } = await apiClient.get(
		`/v1/rooms/${roomId}/default-camera`,
	);
	return data;
};

export const setDefaultCamera = async (
	roomId: string,
	cameraKey: string | null,
): Promise<DefaultCameraResponse> => {
	const { data } = await apiClient.put(
		`/v1/rooms/${roomId}/default-camera`,
		{ default_camera: cameraKey },
	);
	return data;
};

// ==================== Selections API ====================

export interface SelectionGroupResponse {
	group: Record<string, number[]>;
}

export interface SelectionGroupsListResponse {
	items: Record<string, Record<string, number[]>>;
}

export const updateSelection = async (
	roomId: string,
	geometry: string,
	indices: number[],
): Promise<{ status: string }> => {
	const { data } = await apiClient.put(
		`/v1/rooms/${roomId}/geometries/${geometry}/selection`,
		{
			indices,
		},
	);
	return data;
};

export const listSelectionGroups = async (
	roomId: string,
): Promise<SelectionGroupsListResponse> => {
	const { data } = await apiClient.get(
		`/v1/rooms/${roomId}/selection-groups`,
	);
	return data;
};

export const createUpdateSelectionGroup = async (
	roomId: string,
	groupName: string,
	groupData: Record<string, number[]>,
): Promise<{ status: string }> => {
	const { data } = await apiClient.put(
		`/v1/rooms/${roomId}/selection-groups/${groupName}`,
		{ selections: groupData },
	);
	return data;
};

export const deleteSelectionGroup = async (
	roomId: string,
	groupName: string,
): Promise<void> => {
	await apiClient.delete(
		`/v1/rooms/${roomId}/selection-groups/${groupName}`,
	);
};

// ==================== Bookmarks API ====================

export interface BookmarksResponse {
	items: Record<number, string>;
}

export const getAllBookmarks = async (
	roomId: string,
): Promise<BookmarksResponse> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}/bookmarks`);
	// Convert string keys to number keys (JSON serializes number keys as strings)
	const bookmarks: Record<number, string> = {};
	for (const [key, value] of Object.entries(data.items || {})) {
		bookmarks[Number.parseInt(key, 10)] = value as string;
	}
	return { items: bookmarks };
};

export const setBookmark = async (
	roomId: string,
	index: number,
	label: string,
): Promise<{ status: string }> => {
	const { data } = await apiClient.put(
		`/v1/rooms/${roomId}/bookmarks/${index}`,
		{ label },
	);
	return data;
};

export const deleteBookmark = async (
	roomId: string,
	index: number,
): Promise<{ status: string }> => {
	const { data } = await apiClient.delete(
		`/v1/rooms/${roomId}/bookmarks/${index}`,
	);
	return data;
};

// ==================== Room API ====================

export interface CreateRoomRequest {
	room_id: string;
	description?: string;
	copy_from?: string; // Room ID, or @-prefixed preset (@empty, @none)
}

export interface CreateRoomResponse {
	status: string;
	room_id: string;
	frame_count: number;
	created: boolean;
}

export const createRoom = async (
	request: CreateRoomRequest,
	signal?: AbortSignal,
): Promise<CreateRoomResponse> => {
	const { data } = await apiClient.post("/v1/rooms", request, { signal });
	return data;
};

export interface RoomInfo {
	id: string;
	description: string | null;
	frame_count: number;
	locked: boolean;
	hidden: boolean;
}

export const getRoomInfo = async (
	roomId: string,
	signal?: AbortSignal,
): Promise<RoomInfo> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}`, { signal });
	return data;
};

export const getFrameSelection = async (
	roomId: string,
	signal?: AbortSignal,
): Promise<{ frame_selection: number[] | null }> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}/frame-selection`, {
		signal,
	});
	return data;
};

export const updateFrameSelection = async (
	roomId: string,
	indices: number[],
): Promise<{ success: boolean }> => {
	const { data } = await apiClient.put(`/v1/rooms/${roomId}/frame-selection`, {
		indices,
	});
	return data;
};

export const getCurrentStep = async (
	roomId: string,
	signal?: AbortSignal,
): Promise<{ step: number }> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}/step`, { signal });
	return data;
};

// ==================== Jobs API ====================

export interface PaginatedResponse<T> {
	items: T[];
	total: number;
	limit: number;
	offset: number;
}

export interface JobSummary {
	full_name: string;
	category: string;
	name: string;
	workers: string[];
}

export interface JobResponse {
	id: string;
	full_name: string;
	category: string;
	name: string;
	schema: Record<string, unknown>;
	workers: string[];
	room_id: string;
}

export const listJobs = async (
	roomId: string,
): Promise<PaginatedResponse<JobSummary>> => {
	const { data } = await apiClient.get(`/v1/joblib/rooms/${roomId}/jobs`);
	return data;
};

export const getJobSchema = async (
	roomId: string,
	jobName: string,
): Promise<JobResponse> => {
	const { data } = await apiClient.get(
		`/v1/joblib/rooms/${roomId}/jobs/${encodeURIComponent(jobName)}`,
	);
	return data;
};

export interface SubmitTaskResponse {
	id: string;
	job_name: string;
	room_id: string;
	status: TaskStatus;
	created_at: string;
	queue_position: number | null;
}

export const submitTask = async (
	roomId: string,
	jobName: string,
	payload: Record<string, unknown>,
): Promise<SubmitTaskResponse> => {
	const { data } = await apiClient.post(
		`/v1/joblib/rooms/${roomId}/tasks/${encodeURIComponent(jobName)}`,
		{ payload },
	);
	return data;
};

// ==================== Settings API ====================

export interface SettingsResponse {
	schema: any;
	data: Record<string, any>;
}

export const getSettings = async (
	roomId: string,
	sessionId: string,
): Promise<SettingsResponse> => {
	const { data } = await apiClient.get(
		`/v1/rooms/${roomId}/sessions/${sessionId}/settings`,
	);
	return data;
};

export const updateSettings = async (
	roomId: string,
	sessionId: string,
	settingsData: Record<string, any>,
): Promise<{ status: string }> => {
	const { data } = await apiClient.put(
		`/v1/rooms/${roomId}/sessions/${sessionId}/settings`,
		settingsData,
	);
	return data;
};

// ==================== Tasks API ====================

export interface Task {
	id: string;
	job_name: string;
	room_id: string;
	status: TaskStatus;
	created_at: string;
	started_at?: string;
	completed_at?: string;
	queue_position?: number | null;
	error?: string;
	result?: any;
	worker_id?: string | null;
	payload?: Record<string, unknown>;
}

export const listTasksForJob = async (
	roomId: string,
	jobName: string,
	options?: {
		status?: TaskStatus;
		limit?: number;
		offset?: number;
	},
): Promise<PaginatedResponse<Task>> => {
	const params = new URLSearchParams();
	if (options?.status) params.append("status", options.status);
	if (options?.limit !== undefined)
		params.append("limit", options.limit.toString());
	if (options?.offset !== undefined)
		params.append("offset", options.offset.toString());
	const queryString = params.toString();
	const url = `/v1/joblib/rooms/${roomId}/jobs/${encodeURIComponent(jobName)}/tasks${queryString ? `?${queryString}` : ""}`;
	const { data } = await apiClient.get(url);
	return data;
};

export type TypedArray =
	| Float32Array
	| Float64Array
	| Int8Array
	| Int16Array
	| Int32Array
	| BigInt64Array
	| Uint8Array
	| Uint16Array
	| Uint32Array
	| BigUint64Array;

export interface FrameResponse {
	[key: string]: TypedArray;
}

// ==================== Trajectory/Frames API ====================

export const getFrames = async (
	roomId: string,
	frameIndex: number,
	keys: string[],
	signal?: AbortSignal,
): Promise<FrameResponse | null> => {
	const params = new URLSearchParams();
	params.append("indices", frameIndex.toString());
	// because frameIndex is an integer, in decodeTypedData we use [0][key] to access the data
	if (keys && keys.length > 0) {
		params.append("keys", keys.join(","));
	}

	try {
		const response = await apiClient.get(
			`/v1/rooms/${roomId}/frames?${params.toString()}`,
			{
				signal,
				responseType: "arraybuffer",
			},
		);

		const result: Record<string, any> = {};
		for (const key of keys) {
			const decoded = decodeTypedData(response.data, key);
			if (decoded !== undefined) {
				result[key] = decoded;
			}
		}
		return Object.keys(result).length > 0 ? (result as FrameResponse) : null;
	} catch (error: any) {
		if (error?.response?.status === 404) {
			const retryAfter = error.response.headers?.["retry-after"];
			if (retryAfter) {
				// Provider is still computing — throw so fetchWithRetry retries
				throw error;
			}
			// Permanent 404 — frame genuinely doesn't exist
			return null;
		}
		throw error;
	}
};

export interface PartialFrameUpdateResponse {
	success: boolean;
	frame_id: number;
	updated_keys: string[];
}

/**
 * Partially update a frame by merging new data with existing frame data.
 * This allows updating specific keys (e.g., arrays.positions) without sending the entire frame.
 *
 * @param roomId - Room ID
 * @param frameId - Frame index to update
 * @param updates - Dictionary of key-value pairs to update (e.g., {"arrays.positions": Float32Array})
 * @param shapes - Optional map of key paths to shapes for TypedArrays
 * @returns Response with success status and updated keys
 */
export const partialUpdateFrame = async (
	roomId: string,
	frameId: number,
	updates: Record<string, any>,
	shapes?: Map<string, number[]>,
): Promise<PartialFrameUpdateResponse> => {
	const encoded = packBinary(updates, shapes);

	// Use Blob to ensure axios sends exactly the right number of bytes
	// Without this, axios may send extra buffer data beyond the actual content
	// Create a new Uint8Array to get a clean ArrayBuffer (satisfies TypeScript's BlobPart type)
	const blob = new Blob([new Uint8Array(encoded)], {
		type: "application/msgpack",
	});

	const { data } = await apiClient.patch(
		`/v1/rooms/${roomId}/frames/${frameId}`,
		blob,
		{
			headers: {
				"Content-Type": "application/msgpack",
			},
		},
	);

	return data;
};

// ==================== Chat API ====================

export const getChatMessages = async (
	roomId: string,
	limit: number,
	before?: number,
): Promise<ChatMessagesResponse> => {
	const params = new URLSearchParams();
	params.append("limit", limit.toString());
	if (before !== undefined) {
		params.append("before", before.toString());
	}
	const { data } = await apiClient.get(
		`/v1/rooms/${roomId}/chat/messages?${params}`,
	);
	return data;
};

export const createChatMessage = async (
	roomId: string,
	content: string,
): Promise<ChatMessage> => {
	const { data } = await apiClient.post(`/v1/rooms/${roomId}/chat/messages`, {
		content,
	});
	return data;
};

export const editChatMessage = async (
	roomId: string,
	messageId: number,
	content: string,
): Promise<ChatMessage> => {
	const { data } = await apiClient.patch(
		`/v1/rooms/${roomId}/chat/messages/${messageId}`,
		{ content },
	);
	return data;
};

// ==================== Room Management API ====================

export interface Room {
	id: string;
	description?: string | null;
	frame_count: number;
	locked: boolean;
	is_default?: boolean;
	metadata?: Record<string, string>;
}

export interface RoomDetail {
	id: string;
	description?: string | null;
	frame_count: number;
	locked: boolean;
	is_default?: boolean;
}

export interface RoomUpdateRequest {
	description?: string | null;
	locked?: boolean;
}

export interface DefaultRoomResponse {
	room_id: string | null;
}

export interface EditLockResponse {
	locked: boolean;
	lock_token?: string | null;
	user_id?: string | null;
	sid?: string | null;
	msg?: string | null;
	acquired_at?: number | null;
	ttl?: number | null;
}

export const listRooms = async (search?: string): Promise<Room[]> => {
	const params = search ? `?search=${encodeURIComponent(search)}` : "";
	const { data } = await apiClient.get(`/v1/rooms${params}`);
	return data.items;
};

export const getRoom = async (roomId: string): Promise<RoomDetail> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}`);
	return data;
};

export const updateRoom = async (
	roomId: string,
	updates: RoomUpdateRequest,
): Promise<{ status: string }> => {
	const { data } = await apiClient.patch(`/v1/rooms/${roomId}`, updates);
	return data;
};

export const getDefaultRoom = async (): Promise<DefaultRoomResponse> => {
	const { data } = await apiClient.get("/v1/server-settings/default-room");
	return data;
};

export const setDefaultRoom = async (
	roomId: string | null,
): Promise<DefaultRoomResponse> => {
	if (roomId === null) {
		await apiClient.delete("/v1/server-settings/default-room");
		return { room_id: null };
	}
	const { data } = await apiClient.put("/v1/server-settings/default-room", {
		room_id: roomId,
	});
	return data;
};

export const getEditLockStatus = async (
	roomId: string,
): Promise<EditLockResponse> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}/edit-lock`);
	return data;
};

// ==================== Trajectory Upload API ====================

export const uploadTrajectory = async (
	roomId: string,
	file: File,
): Promise<void> => {
	const formData = new FormData();
	formData.append("file", file);
	await apiClient.post(`/v1/rooms/${roomId}/trajectory`, formData, {
		headers: { "Content-Type": "multipart/form-data" },
	});
};

// ==================== Tools API ====================

export interface ConvertMoleculeToImageRequest {
	type: "smiles" | "inchi";
	data: string;
	dark?: boolean;
}

export interface ConvertMoleculeToImageResponse {
	image: string; // data:image/png;base64,...
	status: string;
}

export const convertMoleculeToImage = async (
	request: ConvertMoleculeToImageRequest,
): Promise<ConvertMoleculeToImageResponse> => {
	const { data } = await apiClient.post("/v1/tools/rdkit-img", request);
	return data;
};

export interface DownloadFramesRequest {
	roomId: string;
	indices?: number[];
	selection?: number[];
	filename?: string;
}

export const downloadFrames = async (
	request: DownloadFramesRequest,
): Promise<void> => {
	// 1. Create a temporary download token (authenticated)
	const { data } = await apiClient.post<{ url: string }>(
		`/v1/rooms/${request.roomId}/trajectory/download-tokens`,
	);

	// 2. Build the download URL with query params appended to the token URL
	const tokenUrl = new URL(data.url);
	if (request.indices !== undefined) {
		tokenUrl.searchParams.set("indices", request.indices.join(","));
	}
	if (request.selection) {
		tokenUrl.searchParams.set("selection", request.selection.join(","));
	}
	if (request.filename) {
		tokenUrl.searchParams.set("filename", request.filename);
	}

	// 3. Navigate — browser handles streaming download natively
	window.location.href = tokenUrl.toString();
};

/**
 * Categorize properties based on metadata shape into per-particle and global properties.
 * Per-particle properties have their first dimension equal to the particle count.
 *
 * @param metadata - Frame metadata containing property keys and metadata
 * @param particleCount - Number of particles in the frame
 * @returns Object containing categorized perParticle and global property arrays
 */
export const categorizeProperties = (
	metadata: FrameMetadata,
	particleCount: number,
): {
	perParticle: import("../types/property-inspector").PropertyInfo[];
	global: import("../types/property-inspector").PropertyInfo[];
} => {
	const perParticle: import("../types/property-inspector").PropertyInfo[] = [];
	const global: import("../types/property-inspector").PropertyInfo[] = [];

	for (const [key, meta] of Object.entries(metadata.metadata)) {
		// Per-particle detection: first dimension equals particle count
		const isPerParticle =
			meta.shape &&
			Array.isArray(meta.shape) &&
			meta.shape.length > 0 &&
			meta.shape[0] === particleCount;

		const [prefix] = key.split(".");

		const info: import("../types/property-inspector").PropertyInfo = {
			key,
			metadata: meta,
			category: isPerParticle ? "per-particle" : "global",
			prefix: prefix || "other",
			enabled: false,
		};

		if (isPerParticle) {
			perParticle.push(info);
		} else {
			global.push(info);
		}
	}

	return { perParticle, global };
};

// ==================== Screenshot API ====================

export const completeScreenshot = async (
	uploadUrl: string,
	file: Blob,
	format: string,
	width: number,
	height: number,
): Promise<void> => {
	const formData = new FormData();
	formData.append("file", file, `screenshot.${format}`);
	formData.append("format", format);
	formData.append("width", width.toString());
	formData.append("height", height.toString());

	await apiClient.patch(uploadUrl, formData, {
		headers: { "Content-Type": "multipart/form-data" },
	});
};

// ==================== Server Management API ====================

/**
 * Shutdown the ZnDraw server gracefully.
 */
export const shutdownServer = async (): Promise<{ success: boolean }> => {
	const { data } = await apiClient.post("/v1/admin/shutdown");
	return data;
};

// ==================== Provider / Filesystem APIs ====================

export interface ProviderInfo {
	id: string;
	room_id: string;
	category: string;
	name: string;
	full_name: string;
	schema: Record<string, unknown>;
	worker_id: string;
	created_at: string;
}

export interface FilesystemFileItem {
	name: string;
	path: string;
	size: number;
	type: "file" | "directory";
}

/**
 * List providers registered for a room.
 * Optionally filter by category (e.g., "filesystem").
 */
export const listProviders = async (
	roomId: string,
	category?: string,
): Promise<ProviderInfo[]> => {
	const { data } = await apiClient.get(
		`/v1/joblib/rooms/${roomId}/providers`,
	);
	const items = data.items as ProviderInfo[];
	return category ? items.filter((p) => p.category === category) : items;
};

/**
 * Read data from a provider.
 * The backend long-polls until the result is ready, so this always returns data.
 */
export const readProvider = async (
	roomId: string,
	providerName: string,
	params?: Record<string, string>,
): Promise<FilesystemFileItem[]> => {
	const queryParams = params ? `?${new URLSearchParams(params)}` : "";
	const response = await apiClient.get(
		`/v1/joblib/rooms/${roomId}/providers/${encodeURIComponent(providerName)}${queryParams}`,
	);
	return response.data;
};

// ==================== Step API ====================

export interface StepResponse {
	step: number;
	total_frames: number;
}

export interface UpdateStepRequest {
	step: number;
}

export interface UpdateStepResponse {
	success: boolean;
	step: number;
}

/**
 * Get current step/frame for a room.
 *
 * @param roomId - Room identifier
 * @returns Promise with current step and total frame count
 */
export const getStep = async (roomId: string): Promise<StepResponse> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}/step`);
	return data;
};

/**
 * Set current step/frame for a room.
 * Requires holding the 'step' lock for the room.
 *
 * @param roomId - Room identifier
 * @param step - Frame index to set
 * @returns Promise with update result
 */
export const updateStep = async (
	roomId: string,
	step: number,
): Promise<UpdateStepResponse> => {
	const { data } = await apiClient.put(`/v1/rooms/${roomId}/step`, { step });
	return data;
};

// ==================== Edit Lock API ====================

/**
 * Acquire or refresh the room edit lock.
 *
 * @param roomId - Room identifier
 * @param msg - Optional message describing the lock purpose
 * @param lockToken - Optional lock token for refresh (omit for new acquisition)
 * @returns Promise with edit lock response
 */
export const acquireEditLock = async (
	roomId: string,
	msg?: string,
	lockToken?: string,
): Promise<EditLockResponse> => {
	try {
		const headers: Record<string, string> = {};
		if (lockToken) {
			headers["Lock-Token"] = lockToken;
		}
		const { data } = await apiClient.put(
			`/v1/rooms/${roomId}/edit-lock`,
			{ msg: msg ?? null },
			{ headers },
		);
		return data;
	} catch (error: any) {
		if (error.response?.status === 423) {
			useAppStore
				.getState()
				.showSnackbar(
					"Cannot perform action - room is locked by another user",
					"warning",
				);
		}
		throw error;
	}
};

/**
 * Release the room edit lock.
 *
 * @param roomId - Room identifier
 * @param lockToken - Optional lock token for authentication
 */
export const releaseEditLock = async (
	roomId: string,
	lockToken?: string,
): Promise<void> => {
	const headers: Record<string, string> = {};
	if (lockToken) {
		headers["Lock-Token"] = lockToken;
	}
	await apiClient.delete(`/v1/rooms/${roomId}/edit-lock`, { headers });
};

// ==================== Session API ====================

/**
 * List the current user's active frontend sessions in a room.
 *
 * @param roomId - Room identifier
 * @returns Promise with array of session IDs
 */
export const listSessions = async (roomId: string): Promise<string[]> => {
	const { data } = await apiClient.get(`/v1/rooms/${roomId}/sessions`);
	return data.items;
};

// ==================== Active Camera API ====================

export const updateActiveCamera = async (
	roomId: string,
	sessionId: string,
	cameraKey: string,
): Promise<void> => {
	await apiClient.put(
		`/v1/rooms/${roomId}/sessions/${sessionId}/active-camera`,
		{ active_camera: cameraKey },
	);
};

// ==================== Admin API ====================

export interface AdminUser {
	id: string;
	email: string;
	is_superuser: boolean;
}

export interface AdminUsersResponse {
	items: AdminUser[];
	total: number;
	offset: number;
	limit: number;
}

export const listAdminUsers = async (): Promise<AdminUsersResponse> => {
	const { data } = await apiClient.get("/v1/admin/users");
	return data;
};

export const updateAdminUser = async (
	userId: string,
	update: { is_superuser?: boolean },
): Promise<AdminUser> => {
	const { data } = await apiClient.patch(`/v1/admin/users/${userId}`, update);
	return data;
};

export const deleteAdminUser = async (userId: string): Promise<void> => {
	await apiClient.delete(`/v1/admin/users/${userId}`);
};
