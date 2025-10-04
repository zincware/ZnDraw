import axios from 'axios';

// Define the types based on your Python code
export interface FigureData {
  type: 'plotly';
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

export const listFigures = async (roomId: string): Promise<FigureListResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/figures`);
  return data;
};

export const getFigure = async (roomId:string, key: string): Promise<FigureResponse> => {
  const { data } = await apiClient.get(`/api/rooms/${roomId}/figures/${key}`);
  return data;
};

export const createFigure = async (
  roomId: string,
  key: string,
  figure: FigureData
): Promise<{ status: string }> => {
  const { data } = await apiClient.post(`/api/rooms/${roomId}/figures`, { key, figure });
  return data;
};

export const deleteFigure = async (
  roomId: string,
  key: string
): Promise<{ status: string }> => {
  const { data } = await apiClient.delete(`/api/rooms/${roomId}/figures/${key}`);
  return data;
};