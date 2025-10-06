import { decode } from "@msgpack/msgpack";
import { getFrames } from "../myapi/client";

const numpyDtypeToTypedArray = {
  float32: Float32Array,
  float64: Float64Array,
  int8: Int8Array,
  int16: Int16Array,
  int32: Int32Array,
  uint8: Uint8Array,
  uint16: Uint16Array,
  uint32: Uint32Array,
};

export const getFrameDataOptions = (
  roomId: string,
  frameIndex: number,
  key: string,
) => {
  return {
    // We are not calling useQuery here, but returning its options object
    // so that `useQueries` can use it.
    queryKey: ["frame", roomId, frameIndex, key],
    queryFn: ({ signal }: { signal: AbortSignal }) => getFrames(roomId, frameIndex, [key], signal),
    staleTime: Infinity,
    gcTime: 1000 * 60 * 5,
    // placeholderData: (previousData: any, previousQuery: any) => previousData,
    enabled: !!roomId && frameIndex !== null && !!key,
    select: (arrayBuffer: ArrayBuffer) => {
      const decodedData = decode(arrayBuffer) as any[]; // Type based on your server response
      const singleFrameData = decodedData[0];

      const keyData = singleFrameData[key];
      console.log(`Data for key "${key}" at frame ${frameIndex}:`, keyData);
      if (!keyData) {
        console.warn(
          `Data for key "${key}" not found in response for frame ${frameIndex}`,
        );
        return null;
      }
      if (!keyData.dtype) {
        console.warn(
          `Data for key "${key}" is missing dtype information for frame ${frameIndex}`,
        );
        return null;
      }
      if (!keyData.data) {
        console.warn(
          `Data for key "${key}" is missing actual data for frame ${frameIndex}`,
        );
        return null;
      }

      const TypedArray = numpyDtypeToTypedArray[keyData.dtype as keyof typeof numpyDtypeToTypedArray];
      if (!TypedArray) {
        throw new Error(`Unsupported dtype: ${keyData.dtype}`);
      }

      const dataArray = new TypedArray(keyData.data.slice().buffer);
      if (!dataArray) {
        throw new Error(
          `Failed to create typed array for dtype: ${keyData.dtype}`,
        );
      }
      // console.log(`Fetched key "${key}" for frame ${frameIndex}:`, dataArray);

      return {
        data: dataArray,
        shape: keyData.shape,
      };
    },
  };
};
