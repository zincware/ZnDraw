import axios from "axios";

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
// ==================== Room API ====================

export interface JoinRoomRequest {
  userId: string;
  template?: string;
}

export interface JoinRoomResponse {
  status: string;
  clientId: string;
  joinToken: string;
  frameCount: number;
  roomId: string;
  template: string;
  selection: number[] | null;
  frame_selection: number[] | null;
  created: boolean;
  "presenter-lock": boolean | string | null;
  step: number | null;
  geometries: Record<string, GeometryData> | null;
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

// ==================== Trajectory/Frames API ====================

export const getFrames = async (
  roomId: string,
  frameIndex: number,
  keys: string[],
  signal?: AbortSignal,
): Promise<ArrayBuffer> => {
  const params = new URLSearchParams();
  params.append("indices", frameIndex.toString());
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
  return data;
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

export const listRooms = async (): Promise<Room[]> => {
  const { data } = await apiClient.get("/api/rooms");
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

