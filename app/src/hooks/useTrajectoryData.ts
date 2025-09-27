import { decode } from '@msgpack/msgpack';

const numpyDtypeToTypedArray = {
    'float32': Float32Array,
    'float64': Float64Array,
    'int8': Int8Array,
    'int16': Int16Array,
    'int32': Int32Array,
    'uint8': Uint8Array,
    'uint16': Uint16Array,
    'uint32': Uint32Array,
};


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

export const getFrameDataOptions = (roomId: string, frameIndex: number, key: string) => {

  return {
    // We are not calling useQuery here, but returning its options object
    // so that `useQueries` can use it.
    queryKey: ['frame', roomId, frameIndex, key],
    queryFn: ({ signal }) => fetchFrameData(roomId, frameIndex, [key], signal),
    staleTime: Infinity,
    gcTime: 1000 * 60 * 5,
    // placeholderData: (previousData: any, previousQuery: any) => previousData,
    enabled: !!roomId && frameIndex !== null && !!key,
    select: (arrayBuffer: ArrayBuffer) => {
      const decodedData = decode(arrayBuffer) as any[]; // Type based on your server response
      const singleFrameData = decodedData[0];
      
      const keyData = singleFrameData[key];
      if (!keyData || !keyData.dtype || !keyData.data) {
        console.warn(`Data for key "${key}" not found in response for frame ${frameIndex}`);
        return null;
      }
      
      const TypedArray = numpyDtypeToTypedArray[keyData.dtype];
      if (!TypedArray) {
        throw new Error(`Unsupported dtype: ${keyData.dtype}`);
      }

      const dataArray = new TypedArray(keyData.data.slice().buffer);
      if (!dataArray) {
        throw new Error(`Failed to create typed array for dtype: ${keyData.dtype}`);
      }

      return {
        data: dataArray,
        shape: keyData.shape,
      };
    },
  };
};
