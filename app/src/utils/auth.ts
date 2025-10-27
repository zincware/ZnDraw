/**
 * JWT Authentication utilities for browser client.
 *
 * AUTO-LOGIN STRATEGY:
 * - On first visit, generate a random UUID username
 * - Call /api/login automatically to get JWT token
 * - Store both username and token in localStorage
 * - This avoids username collisions while maintaining seamless UX
 */

const TOKEN_KEY = 'zndraw_jwt_token';
const CLIENT_ID_KEY = 'zndraw_client_id';
const USERNAME_KEY = 'zndraw_username';

export interface LoginResponse {
  status: string;
  token: string;
  clientId: string;
}

/**
 * Get or generate a unique username for this browser.
 * Generates a UUID-based username on first use.
 */
function getOrCreateUsername(): string {
  let username = localStorage.getItem(USERNAME_KEY);

  if (!username) {
    // Generate UUID-based username: "user-abc123..."
    const uuid = crypto.randomUUID();
    username = `user-${uuid.slice(0, 8)}`;
    localStorage.setItem(USERNAME_KEY, username);
    console.log('Generated new username:', username);
  }

  return username;
}

/**
 * Login and get JWT token from server.
 * If userName is not provided, uses the stored/generated username.
 */
export async function login(userName?: string): Promise<LoginResponse> {
  const finalUserName = userName || getOrCreateUsername();

  const response = await fetch('/api/login', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ userName: finalUserName }),
  });

  if (!response.ok) {
    throw new Error(`Login failed: ${response.statusText}`);
  }

  const data = await response.json();

  // Store token and client ID in localStorage
  localStorage.setItem(TOKEN_KEY, data.token);
  localStorage.setItem(CLIENT_ID_KEY, data.clientId);
  if (userName) {
    // If explicit username provided, update stored username
    localStorage.setItem(USERNAME_KEY, userName);
  }

  console.log('Logged in successfully:', data.clientId);

  return data;
}

/**
 * Ensure user is authenticated. Auto-login if needed.
 * Call this before making authenticated requests.
 */
export async function ensureAuthenticated(): Promise<void> {
  const token = getToken();
  if (!token) {
    await login(); // Auto-login with generated username
  }
}

/**
 * Get stored JWT token.
 */
export function getToken(): string | null {
  return localStorage.getItem(TOKEN_KEY);
}

/**
 * Get stored client ID.
 */
export function getClientId(): string | null {
  return localStorage.getItem(CLIENT_ID_KEY);
}

/**
 * Get stored username.
 */
export function getUsername(): string | null {
  return localStorage.getItem(USERNAME_KEY);
}

/**
 * Clear stored authentication data.
 */
export function logout(): void {
  localStorage.removeItem(TOKEN_KEY);
  localStorage.removeItem(CLIENT_ID_KEY);
  localStorage.removeItem(USERNAME_KEY);
  console.log('Logged out');
}

/**
 * Check if user is authenticated (has token).
 */
export function isAuthenticated(): boolean {
  return getToken() !== null;
}
