/**
 * Utility functions for downloading screenshots.
 */

/**
 * Download a screenshot blob to the client's browser.
 *
 * Creates a temporary download link and triggers a file download.
 *
 * @param blob - Screenshot PNG blob
 * @param filename - Optional filename (defaults to timestamp-based name)
 */
export function downloadScreenshot(blob: Blob, filename?: string): void {
	const finalFilename = filename || `screenshot_${Date.now()}.png`;

	// Create object URL for the blob
	const url = URL.createObjectURL(blob);

	// Create temporary download link
	const link = document.createElement("a");
	link.href = url;
	link.download = finalFilename;

	// Trigger download
	document.body.appendChild(link);
	link.click();

	// Cleanup
	document.body.removeChild(link);
	URL.revokeObjectURL(url);
}
