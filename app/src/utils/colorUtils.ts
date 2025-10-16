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
