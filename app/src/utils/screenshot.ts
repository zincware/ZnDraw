/**
 * Utility functions for capturing canvas screenshots.
 * Supports standard canvas capture and is designed to be extensible
 * for future integration with renderers like three-gpu-pathtracer.
 */

/**
 * Capture screenshot from a canvas element as PNG.
 *
 * @param canvas - The canvas element to capture
 * @returns Promise resolving to a PNG Blob
 */
export function captureCanvasScreenshot(
  canvas: HTMLCanvasElement
): Promise<Blob> {
  return new Promise((resolve, reject) => {
    canvas.toBlob(
      (blob) => {
        if (blob) resolve(blob);
        else reject(new Error("Failed to create blob from canvas"));
      },
      "image/png"
    );
  });
}

/**
 * Upload a screenshot blob to the server.
 *
 * @param roomId - Room identifier
 * @param blob - Screenshot PNG blob
 * @param width - Image width
 * @param height - Image height
 * @returns Promise resolving to upload response
 */
export async function uploadScreenshot(
  roomId: string,
  blob: Blob,
  width: number,
  height: number
): Promise<{
  id: string;
  timestamp: string;
  format: string;
  size: number;
  url: string;
}> {
  const formData = new FormData();
  formData.append("file", blob, `screenshot_${Date.now()}.png`);
  formData.append("format", "png");
  formData.append("width", width.toString());
  formData.append("height", height.toString());

  const response = await fetch(`/api/rooms/${roomId}/screenshots/upload`, {
    method: "POST",
    body: formData,
  });

  if (!response.ok) {
    const error = await response.json().catch(() => ({ error: "Upload failed" }));
    throw new Error(error.error || `Upload failed: ${response.statusText}`);
  }

  return response.json();
}

/**
 * Take a screenshot and upload it to the server.
 *
 * @param canvas - Canvas element to capture
 * @param roomId - Room identifier
 * @returns Promise resolving to upload response
 */
export async function takeAndUploadScreenshot(
  canvas: HTMLCanvasElement,
  roomId: string
): Promise<{
  id: string;
  timestamp: string;
  format: string;
  size: number;
  url: string;
}> {
  // Capture screenshot
  const blob = await captureCanvasScreenshot(canvas);

  // Upload to server
  return uploadScreenshot(roomId, blob, canvas.width, canvas.height);
}
