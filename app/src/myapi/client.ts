import axios from "axios";
import { unpackBinary } from "../utils/msgpack-numpy";
import { getToken } from "../utils/auth";
import { useAppStore } from "../store";

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
    const framesList = unpackBinary(encoded instanceof Uint8Array ? encoded : new Uint8Array(encoded));

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
      console.error(`Key "${key}" not found in frame. Available keys:`, Object.keys(frameDict));
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
  figures: string[];
}

export interface GlobalSettings {
  simgen: {
    enabled: boolean;
  };
}

const apiClient = axios.create({});

// Add interceptor to include JWT token and session ID in all requests
apiClient.interceptors.request.use((config) => {
  const token = getToken();
  if (token) {
    config.headers['Authorization'] = `Bearer ${token}`;
  }

  // Add session ID header if available (optional - only needed for specific endpoints)
  // SessionId is required for:
  // - Endpoints with @requires_lock (frame uploads, geometry operations)
  // - Worker registration
  // Most GET endpoints only need JWT auth
  const sessionId = useAppStore.getState().sessionId;
  if (sessionId) {
    config.headers['X-Session-ID'] = sessionId;
  }

  return config;
});

// --- API Functions ---

export const getServerVersion = async (): Promise<{ version: string }> => {
  const { data } = await apiClient.get('/api/version');
  return data;
};

export const getGlobalSettings = async (): Promise<GlobalSettings> => {
  const { data } = await apiClient.get('/api/config/global-settings');
  return data;
};

export const listFigures = async (
  roomId: string,
): Promise<FigureListResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/figures`);
  return data;
};

export const getFigure = async (
  roomId: string,
  key: string,
): Promise<FigureResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/figures/${key}`);
  return data;
};

export const createFigure = async (
  roomId: string,
  key: string,
  figure: FigureData,
): Promise<{ status: string }> => {
  const { data } = await apiClient.post(`/api/rooms/${roomId}/figures`, {
    key,
    figure,
  });
  return data;
};

export const deleteFigure = async (
  roomId: string,
  key: string,
): Promise<{ status: string }> => {
  const { data } = await apiClient.delete(
    `/api/rooms/${roomId}/figures/${key}`,
  );
  return data;
};

export interface FrameKeysResponse {
  frameId: number;
  keys: string[];
  sourceRoom: string;
}

export interface FrameMetadata {
  frameId: number;
  keys: string[];
  metadata: Record<string, any>;
  sourceRoom: string;
}

export const getFrameKeys = async (roomId: string, frameId: number = 0): Promise<FrameKeysResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/frames/${frameId}/keys`);
  return data;
};

export const getFrameMetadata = async (roomId: string, frameId: number = 0): Promise<FrameMetadata> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/frames/${frameId}/metadata`);
  return data;
};

// ==================== Geometries API ====================

export interface GeometryData {
  type: string;
  data: Record<string, any>;
}

export interface GeometryResponse {
  key: string;
  geometry: GeometryData;
}

export interface GeometryData {
  type: string;
  data: any;
}

export interface GeometryListResponse {
  geometries: Record<string, GeometryData>;
}

export const listGeometries = async (
  roomId: string,
): Promise<GeometryListResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/geometries`);
  return data;
};

export const getGeometry = async (
  roomId: string,
  key: string,
): Promise<GeometryResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/geometries/${key}`);
  return data;
};

export const createGeometry = async (
  roomId: string,
  key: string,
  geometryType: string,
  geometryData: Record<string, any>,
  lockToken?: string,
): Promise<{ status: string }> => {
  const requestBody = {
    key,
    type: geometryType,
    data: geometryData,
  };
  
  // If lockToken is provided (from edit mode), client has lock so perform operation directly
  if (lockToken) {
    const { data } = await apiClient.post(`/api/rooms/${roomId}/geometries`, requestBody);
    return data;
  }

  // Otherwise, auto-acquire lock, perform operation, and release
  return withAutoLock(
    roomId,
    "trajectory:meta",
    async (_autoLockToken) => {
      const { data } = await apiClient.post(`/api/rooms/${roomId}/geometries`, requestBody);
      return data;
    },
    `Creating geometry: ${key}`
  );
};

export const deleteGeometry = async (
  roomId: string,
  key: string,
  lockToken?: string,
): Promise<{ status: string }> => {
  // If lockToken is provided (from edit mode), client has lock so perform operation directly
  if (lockToken) {
    const { data } = await apiClient.delete(`/api/rooms/${roomId}/geometries/${key}`);
    return data;
  }

  // Otherwise, auto-acquire lock, perform operation, and release
  return withAutoLock(
    roomId,
    "trajectory:meta",
    async (_autoLockToken) => {
      // Session ID is automatically added by interceptor, server validates session holds lock
      const { data } = await apiClient.delete(`/api/rooms/${roomId}/geometries/${key}`);
      return data;
    },
    `Deleting geometry: ${key}`
  );
};

export const getGeometrySchemas = async (
  roomId: string,
): Promise<{ schemas: Record<string, any> }> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/geometries/schemas`);
  return data;
};

// ==================== Selections API ====================

export interface SelectionsResponse {
  selections: Record<string, number[]>;
  groups: Record<string, Record<string, number[]>>;
  activeGroup: string | null;
}

export interface SelectionResponse {
  geometry: string;
  selection: number[];
}

export interface SelectionGroupResponse {
  group: Record<string, number[]>;
}

export const getAllSelections = async (
  roomId: string,
): Promise<SelectionsResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/selections`);
  return data;
};

export const getSelection = async (
  roomId: string,
  geometry: string,
): Promise<SelectionResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/selections/${geometry}`);
  return data;
};

export const updateSelection = async (
  roomId: string,
  geometry: string,
  indices: number[],
): Promise<{ success: boolean }> => {
  const { data } = await apiClient.put(`/api/rooms/${roomId}/selections/${geometry}`, {
    indices,
  });
  return data;
};

export const getSelectionGroup = async (
  roomId: string,
  groupName: string,
): Promise<SelectionGroupResponse> => {
  const { data } = await apiClient.get(
    `/api/rooms/${roomId}/selections/groups/${groupName}`,
  );
  return data;
};

export const createUpdateSelectionGroup = async (
  roomId: string,
  groupName: string,
  groupData: Record<string, number[]>,
): Promise<{ success: boolean }> => {
  const { data } = await apiClient.put(
    `/api/rooms/${roomId}/selections/groups/${groupName}`,
    groupData,
  );
  return data;
};

export const deleteSelectionGroup = async (
  roomId: string,
  groupName: string,
): Promise<{ success: boolean }> => {
  const { data } = await apiClient.delete(
    `/api/rooms/${roomId}/selections/groups/${groupName}`,
  );
  return data;
};

export const loadSelectionGroup = async (
  roomId: string,
  groupName: string,
): Promise<{ success: boolean }> => {
  const { data } = await apiClient.post(
    `/api/rooms/${roomId}/selections/groups/${groupName}/load`,
  );
  return data;
};

// ==================== Bookmarks API ====================

export interface BookmarksResponse {
  bookmarks: Record<number, string>;
}

export const getAllBookmarks = async (
  roomId: string,
): Promise<BookmarksResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/bookmarks`);
  return data;
};

export const setBookmark = async (
  roomId: string,
  index: number,
  label: string,
): Promise<{ status: string }> => {
  const { data } = await apiClient.put(
    `/api/rooms/${roomId}/bookmarks/${index}`,
    { label },
  );
  return data;
};

export const deleteBookmark = async (
  roomId: string,
  index: number,
): Promise<{ status: string }> => {
  const { data } = await apiClient.delete(
    `/api/rooms/${roomId}/bookmarks/${index}`,
  );
  return data;
};

// ==================== Room API ====================

export interface JoinRoomRequest {
  template?: string;
  allowCreate?: boolean;
}

export interface JoinRoomResponse {
  status: string;
  userName: string;
  sessionId: string; // Session ID for this browser tab
  roomId: string;
  created: boolean;
}

export const joinRoom = async (
  roomId: string,
  request: JoinRoomRequest,
  signal?: AbortSignal,
): Promise<JoinRoomResponse> => {
  const { data } = await apiClient.post(`/api/rooms/${roomId}/join`, request, {
    signal,
  });
  return data;
};

export interface RoomInfo {
  id: string;
  description: string | null;
  frameCount: number;
  locked: boolean;
  hidden: boolean;
}

export const getRoomInfo = async (
  roomId: string,
  signal?: AbortSignal,
): Promise<RoomInfo> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}`, { signal });
  return data;
};

export const getFrameSelection = async (
  roomId: string,
  signal?: AbortSignal,
): Promise<{ frameSelection: number[] | null }> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/frame-selection`, {
    signal,
  });
  return data;
};

export const getCurrentStep = async (
  roomId: string,
  signal?: AbortSignal,
): Promise<{ step: number }> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/step`, { signal });
  return data;
};

export const getUserRoomSettings = async (
  roomId: string,
  signal?: AbortSignal,
): Promise<{ settings: Record<string, any> }> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/settings`, {
    signal,
  });
  return data;
};

// ==================== Templates API ====================

export interface Template {
  id: string;
  name: string;
  description: string;
}

export const listTemplates = async (): Promise<Template[]> => {
  const { data } = await apiClient.get("/api/templates");
  return data;
};

export const getTemplate = async (templateId: string): Promise<Template> => {
  const { data } = await apiClient.get(`/api/templates/${templateId}`);
  return data;
};

// ==================== Schemas/Extensions API ====================

export interface ExtensionMetadata {
  schema: any;
  provider: "celery" | number;
  queueLength: number;
  idleWorkers: number;
  progressingWorkers: number;
}

export interface SchemasResponse {
  [extensionName: string]: ExtensionMetadata;
}

export const getSchemas = async (
  roomId: string,
  category: string,
): Promise<SchemasResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/schema/${category}`);
  return data;
};

export const getExtensionData = async (
  roomId: string,
  category: string,
  extension: string,
): Promise<{ data: any }> => {
  const { data } = await apiClient.get(
    `/api/rooms/${roomId}/extensions/${category}/${extension}/data`,
  );
  return data;
};

export interface SubmitExtensionRequest {
  roomId: string;
  category: string;
  extension: string;
  data: any;
  isPublic?: boolean; // Whether this is a global/public extension (default: false)
}

export const submitExtension = async (
  request: SubmitExtensionRequest,
): Promise<{ status: string; queuePosition?: number; jobId?: string }> => {
  const { roomId, category, extension, data: extensionData, isPublic = false } = request;
  const scope = isPublic ? 'public' : 'private';
  const { data } = await apiClient.post(
    `/api/rooms/${roomId}/extensions/${scope}/${category}/${extension}/submit`,
    {
      data: extensionData,
    },
  );
  return data;
};

// ==================== Jobs API ====================

export interface Job {
  id: string;
  room: string;
  category: string;
  extension: string;
  data: any;
  user_name: string; // snake_case to match backend
  status: "pending" | "assigned" | "processing" | "completed" | "failed";
  provider: string;
  created_at: string;
  assigned_at?: string;
  started_at?: string;
  completed_at?: string;
  worker_id?: string;
  error?: string;
  result?: any;
  queuePosition?: number;
  public?: string; // backend returns this as string "true"/"false"
}

export interface JobsListResponse {
  jobs: Job[];
}

export const listJobs = async (roomId: string): Promise<JobsListResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/jobs`);

  // Backend returns array directly, wrap it in object
  return { jobs: Array.isArray(data) ? data : [] };
};

export const getJob = async (
  roomId: string,
  jobId: string,
): Promise<Job> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/jobs/${jobId}`);
  return data; // API returns job dict directly, not wrapped in {job: ...}
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

  const { data } = await apiClient.get(
    `/api/rooms/${roomId}/frames?${params.toString()}`,
    {
      signal,
      responseType: "arraybuffer",
    },
  );

  const result: Record<string, any> = {};
  for (const key of keys) {
    const decoded = decodeTypedData(data, key);
    if (decoded !== undefined) {
      result[key] = decoded;
    }
  }
  return Object.keys(result).length > 0 ? (result as FrameResponse) : null;
};

// ==================== Chat API ====================

export interface ChatMessage {
  id: string;
  content: string;
  userName: string;
  createdAt: number;
  updatedAt?: number;
  isEdited?: boolean;
}

export interface ChatMetadata {
  hasMore: boolean;
  totalCount: number;
  oldestTimestamp: number | null;
  newestTimestamp: number | null;
}

export interface ChatMessagesResponse {
  messages: ChatMessage[];
  metadata: ChatMetadata;
}

export const getChatMessages = async (
  roomId: string,
  limit: number = 30,
  before?: number,
  after?: number,
): Promise<ChatMessagesResponse> => {
  const params = new URLSearchParams();
  params.append("limit", limit.toString());
  if (before) params.append("before", before.toString());
  if (after) params.append("after", after.toString());

  const { data } = await apiClient.get(
    `/api/rooms/${roomId}/chat/messages?${params.toString()}`,
  );
  return data;
};

// ==================== Room Management API ====================

export interface LockMetadata {
  msg?: string | null;
  userName?: string | null;
  timestamp?: number | null;
}

export interface Room {
  id: string;
  description?: string | null;
  frameCount: number;
  locked: boolean;
  metadataLocked?: LockMetadata | null;
  hidden: boolean;
  isDefault?: boolean;
  presenterSid?: string | null;
  metadata?: Record<string, string>;
}

export interface RoomDetail {
  id: string;
  description?: string | null;
  frameCount: number;
  locked: boolean;
  metadataLocked?: LockMetadata | null;
  hidden: boolean;
  isDefault?: boolean;
  presenterSid?: string | null;
}

export interface RoomUpdateRequest {
  description?: string | null;
  locked?: boolean;
  hidden?: boolean;
}

export interface DuplicateRoomRequest {
  newRoomId?: string;
  description?: string;
}

export interface DuplicateRoomResponse {
  status: string;
  roomId: string;
  frameCount: number;
}

export interface DefaultRoomResponse {
  roomId: string | null;
}

export interface LockStatus {
  locked: boolean;
  target: string;
  holder?: string;
  metadata?: {
    userName?: string;
    msg?: string;
    timestamp?: number;
    [key: string]: any;
  };
  ttl?: number;
}

export const listRooms = async (search?: string): Promise<Room[]> => {
  const params = search ? `?search=${encodeURIComponent(search)}` : "";
  const { data } = await apiClient.get(`/api/rooms${params}`);
  return data;
};

export const getRoom = async (roomId: string): Promise<RoomDetail> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}`);
  return data;
};

export const updateRoom = async (
  roomId: string,
  updates: RoomUpdateRequest,
): Promise<{ status: string }> => {
  const { data } = await apiClient.patch(`/api/rooms/${roomId}`, updates);
  return data;
};

export const duplicateRoom = async (
  roomId: string,
  request: DuplicateRoomRequest = {},
): Promise<DuplicateRoomResponse> => {
  const { data } = await apiClient.post(
    `/api/rooms/${roomId}/duplicate`,
    request,
  );
  return data;
};

export const getDefaultRoom = async (): Promise<DefaultRoomResponse> => {
  const { data } = await apiClient.get("/api/rooms/default");
  return data;
};

export const setDefaultRoom = async (
  roomId: string | null,
): Promise<{ status: string }> => {
  const { data } = await apiClient.put("/api/rooms/default", { roomId });
  return data;
};

export const getLockStatus = async (
  roomId: string,
  target: string = "trajectory:meta",
): Promise<LockStatus> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/locks/${target}`);
  return data;
};

// ==================== File Browser API ====================

export interface FileItem {
  name: string;
  type: "file" | "directory";
  size: number | null;
  modified: string;
  supported?: boolean;
  format_info?: string;
  alreadyLoaded?: {
    room: string;
    description?: string | null;
  };
}

export interface DirectoryListResponse {
  current_path: string;
  items: FileItem[];
  parent: string | null;
}

export interface LoadFileRequest {
  path: string;
  room?: string;
  start?: number;
  stop?: number;
  step?: number;
  make_default?: boolean;
  force_upload?: boolean;
}

export interface LoadFileQueuedResponse {
  status: "queued";
  room: string;
  message: string;
  task_id: string;
}

export interface LoadFileAlreadyLoadedResponse {
  status: "file_already_loaded";
  existingRoom: string;
  message: string;
  filePath: string;
  options: {
    openExisting: string;
    createNew: string;
    forceUpload: string;
  };
}

export type LoadFileResponse = LoadFileQueuedResponse | LoadFileAlreadyLoadedResponse;

export interface SupportedTypesResponse {
  extensions: string[];
  descriptions: Record<string, string>;
}

export interface FileBrowserConfig {
  enabled: boolean;
}

export const listDirectory = async (
  path?: string,
  search?: string,
): Promise<DirectoryListResponse> => {
  const queryParams = new URLSearchParams();
  if (path) queryParams.append("path", path);
  if (search) queryParams.append("search", search);
  const params = queryParams.toString() ? `?${queryParams.toString()}` : "";
  const { data } = await apiClient.get(`/api/file-browser/list${params}`);
  return data;
};

export const loadFile = async (
  request: LoadFileRequest,
): Promise<LoadFileResponse> => {
  const { data } = await apiClient.post("/api/file-browser/load", request);
  return data;
};

export interface UploadFileRequest {
  file: File;
  room?: string;
  start?: number;
  stop?: number;
  step?: number;
  make_default?: boolean;
}

export interface UploadFileResponse {
  status: string;
  room: string;
  task_id: string;
  message: string;
}

export const uploadFile = async (
  request: UploadFileRequest
): Promise<UploadFileResponse> => {
  const formData = new FormData();
  formData.append('file', request.file);

  if (request.room) formData.append('room', request.room);
  if (request.start !== undefined) formData.append('start', request.start.toString());
  if (request.stop !== undefined) formData.append('stop', request.stop.toString());
  if (request.step !== undefined) formData.append('step', request.step.toString());
  if (request.make_default) formData.append('make_default', 'true');

  const { data } = await apiClient.post('/api/file-browser/upload', formData, {
    headers: { 'Content-Type': 'multipart/form-data' }
  });

  return data;
};

export const getSupportedTypes = async (): Promise<SupportedTypesResponse> => {
  const { data } = await apiClient.get("/api/file-browser/supported-types");
  return data;
};

export const getFileBrowserConfig = async (): Promise<FileBrowserConfig> => {
  try {
    // Try to list root directory to check if feature is enabled
    await apiClient.get("/api/file-browser/list");
    return { enabled: true };
  } catch (error) {
    return { enabled: false };
  }
};

export interface CreateRoomFromFileRequest {
  sourceRoom: string;
  newRoom?: string;
  description?: string;
}

export interface CreateRoomFromFileResponse {
  status: string;
  roomId: string;
  sourceRoom: string;
  frameCount: number;
  message: string;
}

export const createRoomFromFile = async (
  request: CreateRoomFromFileRequest,
): Promise<CreateRoomFromFileResponse> => {
  const { data } = await apiClient.post("/api/file-browser/create-room-from-file", request);
  return data;
};

// ==================== Tools API ====================

export interface ConvertMoleculeToImageRequest {
  type: "smiles" | "inchi";
  data: string;
}

export interface ConvertMoleculeToImageResponse {
  image: string; // data:image/png;base64,...
  status: string;
}

export const convertMoleculeToImage = async (
  request: ConvertMoleculeToImageRequest,
): Promise<ConvertMoleculeToImageResponse> => {
  const { data } = await apiClient.post("/api/tools/rdkit-img", request);
  return data;
};

export interface DownloadFramesRequest {
  roomId: string;
  indices?: number[];
  selection?: number[];
  filename?: string;
}

export const downloadFrames = (request: DownloadFramesRequest): void => {
  const params = new URLSearchParams();

  // If indices provided, send them; otherwise no parameters = all frames
  if (request.indices !== undefined) {
    params.append('indices', request.indices.join(','));
  }

  if (request.selection) {
    params.append('selection', request.selection.join(','));
  }

  if (request.filename) {
    params.append('filename', request.filename);
  }

  // Trigger browser download (ExtendedXYZ format)
  const url = `/api/rooms/${request.roomId}/download?${params.toString()}`;
  window.location.href = url;
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
  particleCount: number
): { perParticle: import("../types/property-inspector").PropertyInfo[]; global: import("../types/property-inspector").PropertyInfo[] } => {
  const perParticle: import("../types/property-inspector").PropertyInfo[] = [];
  const global: import("../types/property-inspector").PropertyInfo[] = [];

  metadata.keys.forEach((key) => {
    const meta = metadata.metadata[key];
    if (!meta) return;

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
  });

  return { perParticle, global };
};

// ==================== Server Management API ====================

/**
 * Shutdown the ZnDraw server gracefully.
 */
export const shutdownServer = async (): Promise<{ success: boolean }> => {
  const { data } = await apiClient.post("/api/shutdown");
  return data;
};

// --- Filesystem APIs ---

export interface FilesystemInfo {
  name: string;
  fsType: string;
  public: boolean;
  sessionId: string;
}

export interface FilesystemFileItem {
  name: string;
  path: string;
  size: number;
  type: string;
  modified?: number;
}

export interface ListFilesystemFilesResponse {
  files: FilesystemFileItem[];
}

export interface LoadFilesystemFileRequest {
  path: string;
  targetRoom?: string;
  batchSize?: number;
  start?: number;
  stop?: number;
  step?: number;
}

export interface LoadFilesystemFileResponse {
  success: boolean;
  frameCount: number;
}

export const listFilesystems = async (roomId: string): Promise<FilesystemInfo[]> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/filesystems`);
  return data;
};

export const listFilesystemFiles = async (
  roomId: string,
  fsName: string,
  path?: string,
  recursive?: boolean,
  filterExtensions?: string[]
): Promise<ListFilesystemFilesResponse> => {
  const params = new URLSearchParams();
  if (path) params.append("path", path);
  if (recursive) params.append("recursive", "true");
  if (filterExtensions && filterExtensions.length > 0) {
    params.append("filterExtensions", filterExtensions.join(","));
  }

  const { data } = await apiClient.get(
    `/api/rooms/${roomId}/filesystems/${fsName}/list?${params}`
  );
  return data;
};

export const loadFilesystemFile = async (
  roomId: string,
  fsName: string,
  request: LoadFilesystemFileRequest
): Promise<LoadFilesystemFileResponse> => {
  const { data } = await apiClient.post(
    `/api/rooms/${roomId}/filesystems/${fsName}/load`,
    request
  );
  return data;
};

// Global filesystem operations
export const listGlobalFilesystemFiles = async (
  fsName: string,
  path?: string,
  recursive?: boolean,
  filterExtensions?: string[]
): Promise<ListFilesystemFilesResponse> => {
  const params = new URLSearchParams();
  if (path) params.append("path", path);
  if (recursive) params.append("recursive", "true");
  if (filterExtensions && filterExtensions.length > 0) {
    params.append("filterExtensions", filterExtensions.join(","));
  }

  const { data } = await apiClient.get(
    `/api/filesystems/${fsName}/list?${params}`
  );
  return data;
};

export const loadGlobalFilesystemFile = async (
  fsName: string,
  request: LoadFilesystemFileRequest
): Promise<LoadFilesystemFileResponse> => {
  const { data } = await apiClient.post(
    `/api/filesystems/${fsName}/load`,
    request
  );
  return data;
};

// ==================== Step API ====================

export interface StepResponse {
  step: number;
  totalFrames: number;
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
  const { data } = await apiClient.get(`/api/rooms/${roomId}/step`);
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
  step: number
): Promise<UpdateStepResponse> => {
  const { data } = await apiClient.put(`/api/rooms/${roomId}/step`, { step });
  return data;
};

// ==================== Lock API ====================

export interface LockAcquireResponse {
  success: boolean;
  lockToken?: string;
  ttl?: number;
  refreshInterval?: number;
  error?: string;
}

export interface LockRefreshResponse {
  success: boolean;
  error?: string;
}

export interface LockReleaseResponse {
  success: boolean;
  error?: string;
}

/**
 * Execute an operation with automatic lock acquire/release.
 *
 * This helper acquires a lock, executes the operation with the lock token,
 * and ensures the lock is released even if the operation fails.
 * Shows a snackbar notification if lock acquisition fails.
 *
 * @param roomId - Room identifier
 * @param target - Lock target (e.g., "trajectory:meta")
 * @param operation - Async function that receives the lock token and performs the operation
 * @param msg - Optional message describing the lock purpose
 * @returns Promise with the operation result
 * @throws Error if lock acquisition fails or operation fails
 */
export const withAutoLock = async <T>(
  roomId: string,
  target: string,
  operation: (lockToken: string) => Promise<T>,
  msg?: string
): Promise<T> => {
  // 1. Acquire lock (snackbar shown by acquireLock on failure)
  const acquireResponse = await acquireLock(roomId, target, msg);

  if (!acquireResponse.success || !acquireResponse.lockToken) {
    const errorMsg = acquireResponse.error || "Failed to acquire lock";
    throw new Error(errorMsg);
  }

  const lockToken = acquireResponse.lockToken;

  try {
    // 2. Execute operation with lock token
    const result = await operation(lockToken);

    // 3. Release lock on success
    await releaseLock(roomId, target, lockToken);

    return result;
  } catch (error) {
    // 3. Release lock on failure
    try {
      await releaseLock(roomId, target, lockToken);
    } catch (releaseError) {
      console.error("Failed to release lock after error:", releaseError);
    }

    throw error;
  }
};

/**
 * Acquire a lock for a specific target in a room.
 *
 * @param roomId - Room identifier
 * @param target - Lock target (e.g., "trajectory:meta")
 * @param msg - Optional message describing the lock purpose
 * @returns Promise with lock acquisition result including server-provided TTL and refresh interval
 */
export const acquireLock = async (
  roomId: string,
  target: string,
  msg?: string
): Promise<LockAcquireResponse> => {
  const payload: { msg?: string } = {};
  if (msg) {
    payload.msg = msg;
  }

  try {
    const { data } = await apiClient.post(
      `/api/rooms/${roomId}/locks/${target}/acquire`,
      payload
    );
    return data;
  } catch (error: any) {
    // Show snackbar for lock acquisition failure (423 or other errors)
    useAppStore.getState().showSnackbar(
      "Cannot perform action - room is locked by another user",
      "warning"
    );
    // Re-throw the error so callers know it failed
    throw error;
  }
};

/**
 * Refresh a lock to extend its TTL and optionally update the message.
 *
 * @param roomId - Room identifier
 * @param target - Lock target (e.g., "trajectory:meta")
 * @param lockToken - Lock token from acquire response
 * @param msg - Optional updated message (if provided, updates the lock message)
 * @returns Promise with refresh result
 */
export const refreshLock = async (
  roomId: string,
  target: string,
  lockToken: string,
  msg?: string
): Promise<LockRefreshResponse> => {
  const payload: { lockToken: string; msg?: string } = { lockToken };
  if (msg) {
    payload.msg = msg;
  }

  const { data } = await apiClient.post(
    `/api/rooms/${roomId}/locks/${target}/refresh`,
    payload
  );
  return data;
};

/**
 * Release a lock.
 *
 * @param roomId - Room identifier
 * @param target - Lock target (e.g., "trajectory:meta")
 * @param lockToken - Lock token from acquire response
 * @returns Promise with release result
 */
export const releaseLock = async (
  roomId: string,
  target: string,
  lockToken: string
): Promise<LockReleaseResponse> => {
  const { data } = await apiClient.post(
    `/api/rooms/${roomId}/locks/${target}/release`,
    { lockToken }
  );
  return data;
};

