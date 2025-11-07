/**
 * TypeScript implementation of msgpack-numpy decoder.
 *
 * This is a rewrite of https://github.com/emerald-geomodelling/msgpack-numpy-js
 * with the following improvements:
 * - TypeScript for type safety
 * - Returns multi-dimensional arrays in row-major (C-order) format
 * - Avoids double loops and transposition overhead
 * - Simplified code without external dependencies (underscore)
 */

import { decode as msgpackDecode } from "@msgpack/msgpack";

// Detect system endianness
const isLittleEndian = new Uint8Array(new Uint32Array([0x11223344]).buffer)[0] === 0x44;

// Type mapping from numpy dtype to TypedArray constructor
type TypedArrayConstructor =
  | Uint8ArrayConstructor
  | Int8ArrayConstructor
  | Uint16ArrayConstructor
  | Int16ArrayConstructor
  | Uint32ArrayConstructor
  | Int32ArrayConstructor
  | BigUint64ArrayConstructor
  | BigInt64ArrayConstructor
  | Float32ArrayConstructor
  | Float64ArrayConstructor;

const dtypeToConstructor: Record<string, TypedArrayConstructor> = {
  U8: Uint8Array,
  b: Int8Array,
  B: Uint8Array,
  i1: Int8Array,
  i2: Int16Array,
  i4: Int32Array,
  i8: BigInt64Array,
  I1: Uint8Array,
  I2: Uint16Array,
  I4: Uint32Array,
  I8: BigUint64Array,
  u1: Uint8Array,
  u2: Uint16Array,
  u4: Uint32Array,
  u8: BigUint64Array,
  f4: Float32Array,
  f8: Float64Array,
};

/**
 * Check byte order and swap if needed.
 */
function checkByteOrder(typecode: string, data: Uint8Array): Uint8Array {
  if (typecode === "|") return data; // No byte order
  if (typecode === "=") return data; // Native byte order (assume matches)
  if (typecode === "<" && isLittleEndian) return data;
  if (typecode === ">" && !isLittleEndian) return data;

  throw new Error("Endianness swapping not implemented");
}

/**
 * Reshape a flat TypedArray into multi-dimensional format.
 *
 * msgpack-numpy stores data in C-order (row-major) by default, so we just
 * need to return the flat array as-is for most cases.
 *
 * @param data - Flat typed array with all the data
 * @param shape - Shape tuple [dim0, dim1, ..., dimN]
 * @returns Reshaped array (for 1D/2D, returns flat array; for 3D+, returns nested structure)
 */
function reshapeRowMajor(data: any, shape: number[]): any {
  if (shape.length === 0) {
    // Scalar
    return data[0];
  }

  if (shape.length === 1) {
    // 1D array - return as-is
    return data;
  }

  if (shape.length === 2) {
    // 2D array - most common case for positions, forces, etc.
    // msgpack-numpy stores in C-order (row-major) by default
    // Just verify the data length matches and return as-is

    const [numRows, numCols] = shape;
    const totalElements = numRows * numCols;

    // Verify data length matches shape
    if (data.length !== totalElements) {
      console.warn(`Data length ${data.length} doesn't match shape [${shape.join(", ")}] = ${totalElements}`);
      return data;
    }

    // Data is already in row-major format, return as-is
    return data;
  }

  // For 3D+ arrays, recursively split dimensions
  // Note: This is a simplified implementation. Full n-dimensional transposition
  // would require more complex logic based on memory layout.
  console.warn(`Multi-dimensional arrays (${shape.length}D) not fully optimized, returning as-is`);
  return splitDimensions(data, shape);
}

/**
 * Split flat array into nested array structure (for 3D+ arrays).
 * Falls back to basic dimension splitting for higher-dimensional arrays.
 */
function splitDimensions(data: any, shape: number[]): any {
  if (shape.length === 1) {
    return data;
  }

  const [firstDim, ...restShape] = shape;
  const restSize = restShape.reduce((a, b) => a * b, 1);
  const result: any[] = [];

  for (let i = 0; i < firstDim; i++) {
    const start = i * restSize;
    const end = start + restSize;
    const slice = data.slice(start, end);
    result.push(splitDimensions(slice, restShape));
  }

  return result;
}

/**
 * Unpack msgpack-numpy encoded data.
 *
 * @param value - Decoded msgpack value (may contain numpy arrays or nested msgpack)
 * @returns Unpacked value with TypedArrays in row-major format
 */
function unpackNumpy(value: any): any {
  if (value === null || value === undefined) {
    return value;
  }

  // Check if this is a Uint8Array that might be msgpack-encoded
  // (nested msgpack from the backend)
  if (value instanceof Uint8Array) {
    try {
      // Try to decode as msgpack with numpy support
      const decoded = msgpackDecode(value, {
        mapKeyConverter: (key: unknown) => {
          if (key instanceof Uint8Array) {
            return new TextDecoder().decode(key);
          }
          return key as string | number;
        },
      });
      // Recursively unpack the decoded value
      return unpackNumpy(decoded);
    } catch (err) {
      // If decoding fails, return as-is (might be raw binary data)
      return value;
    }
  }

  if (Array.isArray(value)) {
    return value.map(unpackNumpy);
  }

  if (typeof value === "object") {
    // Check if this is a numpy array encoded as msgpack extension
    // After our mapKeyConverter, keys should be strings
    const isNd = value["nd"];
    const typeStr = value["type"];
    const dataBytes = value["data"];
    const shape = value["shape"];
    const kind = value["kind"];

    if (isNd !== undefined && dataBytes !== undefined) {
      // Handle object dtype (e.g., strings)
      // For object dtype, type is [[b'', b'|O']] and kind is b'O'
      if (kind === "O" || (Array.isArray(typeStr) && typeStr.length > 0)) {
        // Object dtype - data might be pickle-encoded (not supported in JS)
        // msgpack-numpy uses Python pickle for object arrays, which we can't decode in JS
        // The backend should send these as plain lists instead for JS compatibility

        // Check if this looks like pickle data (starts with 0x80 for pickle protocol 2+)
        if (dataBytes instanceof Uint8Array && dataBytes[0] === 0x80) {
          console.warn(
            "Object dtype array uses pickle encoding which is not supported in JavaScript. " +
            "Backend should send string arrays as plain lists instead of numpy object arrays."
          );
          // Return undefined to signal that this data cannot be decoded
          return undefined;
        }

        // Try to decode as msgpack array (for non-pickle object arrays)
        try {
          const decoded = msgpackDecode(dataBytes, {
            mapKeyConverter: (key: unknown) => {
              if (key instanceof Uint8Array) {
                return new TextDecoder().decode(key);
              }
              return key as string | number;
            },
          });
          // Should be an array of strings or other objects
          return decoded;
        } catch (err) {
          console.warn("Failed to decode object dtype array:", err);
          return undefined;
        }
      }

      // Handle unicode string dtype (e.g., <U7)
      if (typeof typeStr === "string" && typeStr[1] === "U") {
        // Unicode string type: <U7 means max 7 characters, UTF-32LE encoded
        const numChars = parseInt(typeStr.slice(2), 10);
        const bytesPerChar = 4; // UTF-32 uses 4 bytes per character
        const bytesPerElement = numChars * bytesPerChar;

        // Check endianness
        const endianChar = typeStr[0];
        if (endianChar !== "<" && endianChar !== "=" && endianChar !== "|") {
          console.warn(`Unsupported endianness for unicode string: ${endianChar}`);
          return undefined;
        }

        // Decode UTF-32LE strings
        const numElements = Array.isArray(shape) ? shape.reduce((a, b) => a * b, 1) : 1;
        const strings: string[] = [];

        for (let i = 0; i < numElements; i++) {
          const offset = i * bytesPerElement;
          const stringBytes = dataBytes.slice(offset, offset + bytesPerElement);

          // Decode UTF-32LE to string
          try {
            // Create a Uint32Array view (UTF-32 = 4 bytes per char)
            const uint32View = new Uint32Array(stringBytes.buffer, stringBytes.byteOffset, numChars);
            // Convert code points to string, stopping at null terminator
            let str = "";
            for (let j = 0; j < numChars; j++) {
              const codePoint = uint32View[j];
              if (codePoint === 0) break; // Null terminator
              str += String.fromCodePoint(codePoint);
            }
            strings.push(str);
          } catch (err) {
            console.warn(`Failed to decode unicode string at index ${i}:`, err);
            strings.push("");
          }
        }

        return strings;
      }

      // Handle numeric dtype
      if (typeof typeStr === "string") {
        const checkedData = checkByteOrder(typeStr[0], dataBytes);
        const dtype = typeStr.slice(1); // Remove endianness prefix

        const TypedArrayCtor = dtypeToConstructor[dtype];
        if (!TypedArrayCtor) {
          console.warn(`Unknown numpy dtype: ${dtype}`);
          return value;
        }

        // Create typed array from the data bytes
        // IMPORTANT: TypedArrays require aligned offsets
        // If the buffer offset is not aligned, we need to copy the data
        let typedArray: any;

        const bytesPerElement = TypedArrayCtor.BYTES_PER_ELEMENT;
        const isAligned = checkedData.byteOffset % bytesPerElement === 0;

        if (isAligned && checkedData.byteOffset !== undefined) {
          // Can create a view directly (aligned offset)
          typedArray = new TypedArrayCtor(
            checkedData.buffer,
            checkedData.byteOffset,
            checkedData.byteLength / bytesPerElement
          );
        } else {
          // Need to copy to a new aligned buffer
          const alignedBuffer = new Uint8Array(checkedData).buffer;
          typedArray = new TypedArrayCtor(alignedBuffer);
        }

        // Reshape to proper dimensions in row-major format
        if (shape && Array.isArray(shape)) {
          return reshapeRowMajor(typedArray, shape);
        }

        return typedArray;
      }
    }

    // Not a numpy array, recursively unpack object properties
    const result: Record<string, any> = {};
    for (const key in value) {
      result[key] = unpackNumpy(value[key]);
    }
    return result;
  }

  return value;
}

/**
 * Decode msgpack binary data with numpy array support.
 * Returns multi-dimensional arrays in row-major (C-order) format.
 *
 * @param data - Msgpack-encoded binary data
 * @returns Decoded data with TypedArrays
 */
/**
 * Convert bytes keys to string keys recursively.
 */
function convertBytesKeys(value: any): any {
  if (value === null || value === undefined) {
    return value;
  }

  if (Array.isArray(value)) {
    return value.map(convertBytesKeys);
  }

  if (typeof value === "object") {
    if (value instanceof Uint8Array || ArrayBuffer.isView(value)) {
      // Leave TypedArrays as-is
      return value;
    }

    // Convert object keys from bytes to strings
    const result: Record<string, any> = {};
    for (const key in value) {
      result[key] = convertBytesKeys(value[key]);
    }
    return result;
  }

  return value;
}

export function unpackBinary(data: Uint8Array | ArrayBuffer): any {
  // Decode msgpack with custom mapKeyConverter to handle bytes keys
  // msgpack-numpy uses bytes keys for metadata like b"nd", b"type", b"data", b"shape"
  const decoded = msgpackDecode(data instanceof Uint8Array ? data : new Uint8Array(data), {
    mapKeyConverter: (key: unknown) => {
      // Convert Uint8Array keys to strings
      if (key instanceof Uint8Array) {
        return new TextDecoder().decode(key);
      }
      // Return other keys as-is
      return key as string | number;
    },
  });

  return unpackNumpy(decoded);
}
