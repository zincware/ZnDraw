/**
 * Layout constants for consistent spacing across the application.
 * These ensure proper alignment of fixed UI elements.
 */

export const LAYOUT_CONSTANTS = {
  /** Height of the top AppBar in pixels */
  APPBAR_HEIGHT: 50,

  /** Height of the bottom FrameProgressBar in pixels */
  PROGRESSBAR_HEIGHT: 45,

  /** Width of the primary drawer (icon sidebar) in pixels */
  PRIMARY_DRAWER_WIDTH: 50,

  /** Width of the secondary drawer (modifiers, settings, etc.) in pixels */
  SECONDARY_DRAWER_WIDTH: 600,
} as const;

/**
 * Calculate the available content height between AppBar and ProgressBar
 */
export const getContentHeight = () => {
  return `calc(100vh - ${LAYOUT_CONSTANTS.APPBAR_HEIGHT}px - ${LAYOUT_CONSTANTS.PROGRESSBAR_HEIGHT}px)`;
};
