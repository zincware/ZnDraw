/**
 * Checks if a value is a hex color string (e.g., "#FF5733" or "#abc")
 * @param value The value to check
 * @returns True if the value is a hex color string
 */
export const isHexColor = (value: any): boolean => {
  if (typeof value !== "string") return false;
  // Match #RGB or #RRGGBB format
  return /^#[0-9A-Fa-f]{3}$|^#[0-9A-Fa-f]{6}$/.test(value);
};

/**
 * Checks if a value should be fetched from the server as frame data.
 * Returns false for hex colors, true for data keys.
 * @param value The value to check
 * @returns True if the value should be fetched from server
 */
export const shouldFetchAsFrameData = (value: any): boolean => {
  if (typeof value !== "string") return false;
  return !isHexColor(value);
};

/**
 * Converts a hex color to RGB array in 0-1 range
 * @param hex The hex color string (e.g., "#FF5733" or "#f73")
 * @returns Array of [r, g, b] in 0-1 range, or null if invalid
 */
export const hexToRgb = (hex: string): [number, number, number] | null => {
  if (!isHexColor(hex)) return null;

  const hexValue = hex.replace("#", "");

  let r: number, g: number, b: number;

  if (hexValue.length === 3) {
    // Convert #RGB to #RRGGBB
    r = parseInt(hexValue[0] + hexValue[0], 16);
    g = parseInt(hexValue[1] + hexValue[1], 16);
    b = parseInt(hexValue[2] + hexValue[2], 16);
  } else {
    r = parseInt(hexValue.slice(0, 2), 16);
    g = parseInt(hexValue.slice(2, 4), 16);
    b = parseInt(hexValue.slice(4, 6), 16);
  }

  // Return in 0-1 range for Three.js
  return [r / 255, g / 255, b / 255];
};
