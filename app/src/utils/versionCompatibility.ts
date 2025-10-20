/**
 * Parses a semantic version string into its components.
 * Handles pre-release versions like "0.6.0a4"
 */
export function parseVersion(version: string): { major: number; minor: number; patch: number } | null {
  // Strip any pre-release suffix (e.g., "a4", "b1", "rc2")
  const cleanVersion = version.replace(/[a-zA-Z].*/g, '');
  const parts = cleanVersion.split('.').map(Number);

  if (parts.length < 3 || parts.some(isNaN)) {
    return null;
  }

  return {
    major: parts[0],
    minor: parts[1],
    patch: parts[2]
  };
}

/**
 * Checks version compatibility between client and server.
 *
 * @param clientVersion - The client version string
 * @param serverVersion - The server version string
 * @returns An object with:
 *   - compatible: boolean indicating if versions are compatible
 *   - severity: 'none' | 'warning' | 'error'
 *   - message: string with human-readable message
 */
export function checkVersionCompatibility(
  clientVersion: string,
  serverVersion: string
): { compatible: boolean; severity: 'none' | 'warning' | 'error'; message: string } {
  const client = parseVersion(clientVersion);
  const server = parseVersion(serverVersion);

  if (!client || !server) {
    return {
      compatible: false,
      severity: 'error',
      message: `Invalid version format. Client: ${clientVersion}, Server: ${serverVersion}`
    };
  }

  // Check major version mismatch
  if (client.major !== server.major) {
    return {
      compatible: false,
      severity: 'error',
      message: `Major version mismatch. Client: ${clientVersion}, Server: ${serverVersion}. Connection refused.`
    };
  }

  // Check minor version mismatch
  if (client.minor !== server.minor) {
    return {
      compatible: false,
      severity: 'error',
      message: `Minor version mismatch. Client: ${clientVersion}, Server: ${serverVersion}. Connection refused.`
    };
  }

  // Check patch version mismatch (warning only)
  if (client.patch !== server.patch) {
    return {
      compatible: true,
      severity: 'warning',
      message: `Patch version mismatch. Client: ${clientVersion}, Server: ${serverVersion}. This may cause unexpected behavior.`
    };
  }

  return {
    compatible: true,
    severity: 'none',
    message: 'Version compatible'
  };
}

/**
 * Gets the client version.
 * The frontend version should match the backend ZnDraw version.
 * This is injected at build time from the Python package version.
 */
export function getClientVersion(): string {
  // This environment variable should be set during the build process
  // from the Python package's __version__
  return import.meta.env.VITE_APP_VERSION || '0.6.0a4';
}
