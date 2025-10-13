/**
 * Utility functions for bookmark key conversion.
 * JSON serialization converts integer keys to strings, so we need to convert them back.
 */

/**
 * Convert bookmark string keys from JSON to number keys.
 *
 * @param bookmarks - Bookmarks object with string keys from JSON
 * @returns Bookmarks object with number keys, or null if input is null/undefined
 *
 * @example
 * convertBookmarkKeys({"1": "First", "5": "Middle"}) // {1: "First", 5: "Middle"}
 */
export function convertBookmarkKeys(
  bookmarks: Record<string, string> | null | undefined
): Record<number, string> | null {
  if (!bookmarks) return null;

  const result: Record<number, string> = {};
  for (const [key, value] of Object.entries(bookmarks)) {
    result[parseInt(key, 10)] = value;
  }
  return result;
}
