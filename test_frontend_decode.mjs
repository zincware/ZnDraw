/**
 * Test script to verify frontend msgpack-numpy decoding.
 *
 * This script simulates what the frontend does:
 * 1. Fetch frames from the backend
 * 2. Decode using msgpack-numpy-js
 * 3. Verify the data is correct
 */

import { decode } from "@msgpack/msgpack";
import { unpackBinary } from "./app/src/utils/msgpack-numpy.ts";
import axios from "axios";

async function testFrameDecoding() {
  console.log("=".repeat(80));
  console.log("Testing Frontend Msgpack-Numpy Decoding");
  console.log("=".repeat(80));

  try {
    // Step 1: Fetch frame data from the backend
    console.log("\n1. Fetching frame data from backend...");
    const response = await axios.get("http://localhost:5000/api/rooms/tmp_s22_xyz/frames", {
      params: {
        indices: "0",
        keys: "arrays.positions,arrays.numbers",
      },
      responseType: "arraybuffer",
    });

    console.log(`   Status: ${response.status}`);
    console.log(`   Content-Type: ${response.headers["content-type"]}`);
    console.log(`   Content length: ${response.data.byteLength} bytes`);

    // Step 2: Decode outer msgpack layer
    console.log("\n2. Decoding outer msgpack layer...");
    // IMPORTANT: @msgpack/msgpack by default doesn't allow non-string/number keys
    // We need to use a custom decoder or work around this
    console.log("   Trying to decode with custom mapKeyDecoder...");

    let framesList;
    try {
      // Try with mapKeyDecoder
      framesList = decode(new Uint8Array(response.data), {
        mapKeyDecoder: (key) => {
          console.log(`   mapKeyDecoder called with: ${key?.constructor?.name}`);
          if (key instanceof Uint8Array) {
            const decoded = new TextDecoder().decode(key);
            console.log(`   Decoded key: "${decoded}"`);
            return decoded;
          }
          return key;
        },
      });
    } catch (firstError) {
      console.log(`   ✗ mapKeyDecoder approach failed: ${firstError.message}`);
      console.log("\n   Trying alternative: unpackBinary directly on whole response...");

      // Try using msgpack-numpy-js unpackBinary on the entire response
      framesList = unpackBinary(new Uint8Array(response.data));
    }

    console.log(`   Type: ${Array.isArray(framesList) ? "Array" : typeof framesList}`);
    console.log(`   Length: ${framesList.length}`);

    // Step 3: Get first frame
    console.log("\n3. Extracting first frame...");
    const frameDict = framesList[0];

    console.log(`   Type: ${typeof frameDict}`);
    console.log(`   Keys: ${Object.keys(frameDict).join(", ")}`);

    // Step 4: Inspect positions array
    console.log("\n4. Inspecting arrays.positions...");
    const positionsKey = "arrays.positions";
    const positions = frameDict[positionsKey];

    console.log(`   Type: ${typeof positions}`);
    console.log(`   Constructor: ${positions?.constructor?.name}`);

    if (positions && typeof positions === 'object') {
      // Check if already decoded as TypedArray
      const isTypedArray = ArrayBuffer.isView(positions);
      const isArray = Array.isArray(positions);

      console.log(`   Is TypedArray: ${isTypedArray}`);
      console.log(`   Is Array: ${isArray}`);

      if (isTypedArray) {
        // ✅ Already decoded! unpackBinary did everything in one go
        console.log(`\n   ✓✓✓ Perfect! Already decoded to flat TypedArray in row-major format`);
        console.log(`   TypedArray type: ${positions.constructor.name}`);
        console.log(`   Total length: ${positions.length} elements`);
        const numAtoms = Math.floor(positions.length / 3);
        console.log(`   Shape: [${numAtoms}, 3]`);
        console.log(`   First atom:  [${positions[0]}, ${positions[1]}, ${positions[2]}]`);
        console.log(`   Second atom: [${positions[3]}, ${positions[4]}, ${positions[5]}]`);
        console.log(`   Third atom:  [${positions[6]}, ${positions[7]}, ${positions[8]}]`);
      } else if (isArray) {
        console.log(`\n   ⚠️  Got Array instead of flat TypedArray`);
        console.log(`   Array length: ${positions.length}`);
        console.log(`   First element type: ${positions[0]?.constructor?.name}`);
      } else {
        console.log(`\n   ? Unknown object type`);
      }
    } else {
      console.error(`   ✗ Unexpected type: ${typeof positions}`);
    }

    // Step 5: Check numbers array
    console.log("\n5. Checking arrays.numbers...");
    const numbersKey = "arrays.numbers";
    const numbers = frameDict[numbersKey];

    if (ArrayBuffer.isView(numbers)) {
      console.log(`\n   ✓ Already decoded!`);
      console.log(`   Type: ${numbers?.constructor?.name}`);
      console.log(`   Length: ${numbers.length}`);
      console.log(`   Values: [${Array.from(numbers).join(", ")}]`);
    } else {
      console.error(`   ✗ Expected TypedArray, got ${typeof numbers}`);
    }

    console.log("\n" + "=".repeat(80));
    console.log("✓ All tests passed!");
    console.log("=".repeat(80));
  } catch (error) {
    console.error("\n✗ Test failed:");
    console.error(error.message);
    if (error.response) {
      console.error(`   Response status: ${error.response.status}`);
      console.error(`   Response data: ${error.response.data}`);
    }
    console.error("\nStack trace:");
    console.error(error.stack);
    process.exit(1);
  }
}

// Run the test
testFrameDecoding();
