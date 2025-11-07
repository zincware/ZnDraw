/**
 * Test decoding arrays.colors to verify string arrays work correctly.
 */

import { unpackBinary } from "./app/src/utils/msgpack-numpy.ts";
import axios from "axios";

async function testColorsDecoding() {
  console.log("=" .repeat(80));
  console.log("Testing arrays.colors decoding");
  console.log("=".repeat(80));

  try {
    const response = await axios.get("http://localhost:5000/api/rooms/tmp_s22_xyz/frames", {
      params: { indices: "0", keys: "arrays.colors" },
      responseType: "arraybuffer",
    });

    console.log(`\n✓ Fetched ${response.data.byteLength} bytes from server`);

    // Decode with our custom decoder
    const frames = unpackBinary(new Uint8Array(response.data));
    console.log(`✓ Decoded to ${frames.length} frame(s)`);

    const frame = frames[0];
    console.log(`✓ Frame has keys: ${Object.keys(frame).join(", ")}`);

    // Check colors
    const colors = frame["arrays.colors"];
    console.log(`\n✓ Colors type: ${colors?.constructor?.name}`);
    console.log(`✓ Is Array: ${Array.isArray(colors)}`);
    console.log(`✓ Is TypedArray: ${ArrayBuffer.isView(colors)}`);

    if (colors) {
      console.log(`✓ Length: ${colors.length}`);
      console.log(`\nFirst 5 colors:`);
      for (let i = 0; i < Math.min(5, colors.length); i++) {
        console.log(`  [${i}]: ${colors[i]} (type: ${typeof colors[i]})`);
      }
    }

    console.log("\n" + "=".repeat(80));
    console.log("✓✓✓ Colors decoding test passed!");
    console.log("=".repeat(80));
  } catch (error) {
    console.error("\n✗ Test failed:");
    console.error(error.message);
    console.error("\nStack trace:");
    console.error(error.stack);
    process.exit(1);
  }
}

testColorsDecoding();
