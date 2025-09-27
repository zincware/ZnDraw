import { useQuery } from '@tanstack/react-query';
import { useAppStore } from '../store';
import { decode } from '@msgpack/msgpack';



interface DecodedNumpyArray {
  shape: number[];
  dtype: string;
  data: any; // This can be more specific based on your needs
}

interface FrameResponseData {
  positions?: DecodedNumpyArray;
  colors?: DecodedNumpyArray;
  radii?: DecodedNumpyArray;
}

const fetchFrameData = async (roomId: string, frameIndex: number, keys: string[], signal: AbortSignal) => {
  const response = await fetch(`/api/frames/${roomId}`, {
    signal, // TanStack Query passes an AbortSignal automatically!
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ indices: [frameIndex], keys: keys }),
  });

  if (!response.ok) {
    console.error("Received error response from /api/frames:", response);
    throw new Error(`Request failed with status ${response.status}`);
  }
  return response.arrayBuffer();
};

// The new hook is much simpler.
export const useFrameData = (frameIndex: number, keys: string[]) => {
  const { roomId } = useAppStore();

  return useQuery({
    // 1. The Query Key: Uniquely identifies this data.
    // When frameIndex changes, TanStack Query will automatically fetch new data.
    queryKey: ['frame', roomId, frameIndex, keys],
    
    // 2. The Query Function: The async function that fetches the data.
    queryFn: ({ signal }) => fetchFrameData(roomId, frameIndex, keys, signal),

    // 3. Configuration Options: This is where the magic happens.
    staleTime: Infinity,   // Frames are immutable, so they never become "stale".
    gcTime: 1000 * 60 * 5, // Garbage collect unused frames after 5 minutes.
    enabled: !!roomId && frameIndex !== null, // Only run the query if we have a roomId and frameIndex.
    
    // 4. Select: Transform the data AFTER it's fetched and cached.
    // This is very efficient. The raw ArrayBuffer is cached, and the conversion
    // only runs when a component needs the data.
    select: (arrayBuffer: ArrayBuffer) => {
      const decodedData = decode(arrayBuffer) as FrameResponseData[];
      const singleFrameData = decodedData[0];
      // assume all are float32
      // TODO: get dtype from key.dtype!
      const positions = new Float64Array(singleFrameData.positions.data.slice().buffer);
      const colors = new Float16Array(singleFrameData.colors.data.slice().buffer);
      const radii = new Float32Array(singleFrameData.radii.data.slice().buffer);

      console.log(singleFrameData);

      return {
        positions, // Float32Array of [x1, y1, z1, x2, y2, z2, ...]
        colors,    // Float32Array of [r1, g1, b1, r2, g2, b2, ...]
        radii,     // Float32Array of [r1, r2, r3, ...]
        count: singleFrameData.positions.shape[0], // The number of atoms
      };
    },
  });
};
