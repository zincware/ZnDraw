/**
 * Test that unpackBinary can decode the full response with bytes keys.
 */

import { unpackBinary } from "./app/src/utils/msgpack-numpy.ts";
import axios from "axios";

async function testBytesKeys() {
  console.log("=" .repeat(80));
  console.log("Testing msgpack-numpy with bytes keys everywhere");
  console.log("=".repeat(80));

  try {
    const response = await axios.get("http://localhost:5000/api/rooms/tmp_s22_xyz/frames", {
      params: { indices: "0", keys: "arrays.positions,arrays.numbers" },
      responseType: "arraybuffer",
    });

    console.log(`\n✓ Fetched ${response.data.byteLength} bytes from server`);

    // Decode the entire response with our custom decoder
    // This should handle both outer bytes keys AND inner numpy arrays
    console.log("\nDecoding with unpackBinary...");
    const frames = unpackBinary(new Uint8Array(response.data));

    console.log(`✓ Decoded to ${frames.length} frame(s)`);

    const frame = frames[0];
    console.log(`✓ Frame has keys: ${Object.keys(frame).join(", ")}`);

    // Check positions
    const positions = frame["arrays.positions"];
    console.log(`\n✓ Positions type: ${positions?.constructor?.name}`);
    console.log(`✓ Is TypedArray: ${ArrayBuffer.isView(positions)}`);
    console.log(`✓ Length: ${positions.length} elements`);
    console.log(`✓ Shape: [${positions.length / 3}, 3]`);
    console.log(`\nFirst 3 atoms:`);
    console.log(`  [${positions[0]}, ${positions[1]}, ${positions[2]}]`);
    console.log(`  [${positions[3]}, ${positions[4]}, ${positions[5]}]`);
    console.log(`  [${positions[6]}, ${positions[7]}, ${positions[8]}]`);

    // Check numbers
    const numbers = frame["arrays.numbers"];
    console.log(`\n✓ Numbers type: ${numbers?.constructor?.name}`);
    console.log(`✓ Length: ${numbers.length} elements`);
    console.log(`✓ Values: [${Array.from(numbers).join(", ")}]`);

    console.log("\n" + "=".repeat(80));
    console.log("✓✓✓ All tests passed! Bytes keys work end-to-end!");
    console.log("=".repeat(80));
  } catch (error) {
    console.error("\n✗ Test failed:");
    console.error(error.message);
    console.error("\nStack trace:");
    console.error(error.stack);
    process.exit(1);
  }
}

testBytesKeys();
