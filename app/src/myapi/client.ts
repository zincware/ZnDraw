import axios from "axios";
import { decode } from "@msgpack/msgpack";

const numpyDtypeToTypedArray = {
  float32: Float32Array,
  float64: Float64Array,
  int8: Int8Array,
  int16: Int16Array,
  int32: Int32Array,
  int64: BigInt64Array,
  uint8: Uint8Array,
  uint16: Uint16Array,
  uint32: Uint32Array,
  uint64: BigUint64Array,
};

function decodeTypedData(encoded: any, key: string) {
  if (!encoded) return undefined;
  try {
    const decoded = decode(encoded)[0][key] as { dtype: string; data: any };

    // Handle object dtype (e.g., hex color strings)
    // Object dtype arrays are sent as JSON-encoded strings in msgpack
    if (decoded.dtype === 'object') {
      // The data field contains a JSON string that needs to be parsed
      if (typeof decoded.data === 'string') {
        return JSON.parse(decoded.data);
      }
      // If already parsed (shouldn't happen, but handle gracefully)
      return decoded.data;
    }

    const TypedArrayCtor = numpyDtypeToTypedArray[decoded.dtype as keyof typeof numpyDtypeToTypedArray];
    if (!TypedArrayCtor) throw new Error(`Unsupported dtype: ${decoded.dtype}`);
    return new TypedArrayCtor(decoded.data.slice().buffer);
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

const apiClient = axios.create({});

// --- API Functions ---

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

export interface FrameMetadata {
  frameId: number;
  keys: string[];
  metadata: Record<string, any>;
  sourceRoom: string;
}

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

export interface GeometryListResponse {
  geometries: string[];
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
  clientId: string | null,
  key: string,
  geometryType: string,
  geometryData: Record<string, any>,
): Promise<{ status: string }> => {
  const { data } = await apiClient.post(`/api/rooms/${roomId}/geometries`, {
    key,
    type: geometryType,
    data: geometryData,
    clientId,
  });
  return data;
};

export const deleteGeometry = async (
  roomId: string,
  key: string,
): Promise<{ status: string }> => {
  const { data } = await apiClient.delete(`/api/rooms/${roomId}/geometries/${key}`);
  return data;
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
  userId: string;
  template?: string;
  allowCreate?: boolean;
}

export interface JoinRoomResponse {
  status: string;
  clientId: string;
  joinToken: string;
  frameCount: number;
  roomId: string;
  template: string;
  selections: Record<string, number[]>;
  selectionGroups: Record<string, Record<string, number[]>>;
  activeSelectionGroup: string | null;
  frame_selection: number[] | null;
  created: boolean;
  "presenter-lock": boolean | string | null;
  step: number | null;
  geometries: Record<string, GeometryData> | null;
  geometryDefaults: Record<string, any> | null;
  bookmarks: Record<number, string> | null;
  settings: Record<string, any>;
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
  userId: string,
  category: string,
  extension: string,
): Promise<{ data: any }> => {
  const { data } = await apiClient.get(
    `/api/rooms/${roomId}/extensions/${category}/${extension}/data`,
    {
      params: { userId },
    },
  );
  return data;
};

export interface SubmitExtensionRequest {
  roomId: string;
  userId: string;
  category: string;
  extension: string;
  data: any;
}

export const submitExtension = async (
  request: SubmitExtensionRequest,
): Promise<{ status: string; queuePosition?: number; jobId?: string }> => {
  const { roomId, userId, category, extension, data: extensionData } = request;
  const { data } = await apiClient.post(
    `/api/rooms/${roomId}/extensions/${category}/${extension}/submit`,
    {
      userId,
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
  user_id: string;
  status: "queued" | "running" | "completed" | "failed";
  provider: string;
  created_at: string;
  started_at: string;
  completed_at: string;
  worker_id: string;
  error: string;
  result: any;
}

export interface JobsListResponse {
  jobs: Job[];
}

export const listJobs = async (roomId: string): Promise<JobsListResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/jobs`);
  return data;
};

export interface JobResponse {
  job: Job;
}

export const getJob = async (
  roomId: string,
  jobId: string,
): Promise<JobResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/jobs/${jobId}`);
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
  userId: string;
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

export interface Room {
  id: string;
  description?: string | null;
  frameCount: number;
  locked: boolean;
  metadataLocked?: boolean;
  hidden: boolean;
  isDefault?: boolean;
  metadata?: Record<string, string>;
}

export interface RoomDetail {
  id: string;
  description?: string | null;
  frameCount: number;
  locked: boolean;
  metadataLocked?: boolean;
  hidden: boolean;
  isDefault?: boolean;
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

