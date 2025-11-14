/**
 * JWT Authentication utilities for browser client.
 *
 * Simplified Architecture:
 * - userName is the primary identifier (unique, immutable after registration)
 * - Guests get auto-generated usernames (user-xyz)
 * - Registration allows choosing permanent username
 * - No client_id concept
 */

const TOKEN_KEY = 'zndraw_jwt_token';
const USERNAME_KEY = 'zndraw_username';
const IS_ADMIN_KEY = 'zndraw_is_admin';
const USER_ROLE_KEY = 'zndraw_user_role';

export type UserRole = 'guest' | 'user' | 'admin';

export interface LoginResponse {
  status: string;
  token: string;
  userName: string;
  role: UserRole;
}

export interface UserRoleResponse {
  userName: string;
  role: UserRole;
}

export interface RegisterResponse {
  status: string;
  token: string;
  userName: string;
  role: UserRole;
}

/**
 * Login and get JWT token from server.
 *
 * @param userName - Optional username (if not provided, server generates guest username)
 * @param password - Optional password (if provided, authenticates as registered user)
 */
export async function login(userName?: string, password?: string): Promise<LoginResponse> {
  const body: { userName?: string; password?: string } = {};

  if (userName !== undefined) {
    body.userName = userName;
  }
  if (password !== undefined) {
    body.password = password;
  }

  const response = await fetch('/api/login', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });

  if (!response.ok) {
    const errorData = await response.json().catch(() => ({ error: response.statusText }));
    throw new Error(errorData.error || `Login failed: ${response.statusText}`);
  }

  const data = await response.json();

  // Store token, username, and role in localStorage
  localStorage.setItem(TOKEN_KEY, data.token);
  localStorage.setItem(USERNAME_KEY, data.userName);
  localStorage.setItem(USER_ROLE_KEY, data.role);

  return data;
}

/**
 * Ensure user is authenticated. Auto-login if needed.
 * Call this before making authenticated requests.
 */
export async function ensureAuthenticated(): Promise<void> {
  const token = getToken();
  if (!token) {
    await login(); // Auto-login with server-generated guest username
  }
}

/**
 * Get stored JWT token.
 */
export function getToken(): string | null {
  return localStorage.getItem(TOKEN_KEY);
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
  localStorage.removeItem(USERNAME_KEY);
  localStorage.removeItem(USER_ROLE_KEY);
}

/**
 * Check if user is authenticated (has token).
 */
export function isAuthenticated(): boolean {
  return getToken() !== null;
}

/**
 * Fetch user role from server.
 */
export async function fetchUserRole(): Promise<UserRoleResponse> {
  const token = getToken();
  if (!token) {
    throw new Error('Not authenticated');
  }

  const response = await fetch('/api/user/role', {
    method: 'GET',
    headers: {
      'Authorization': `Bearer ${token}`,
      'Content-Type': 'application/json',
    },
  });

  if (!response.ok) {
    const errorData = await response.json().catch(() => ({ error: response.statusText }));
    throw new Error(errorData.error || `Failed to fetch role: ${response.statusText}`);
  }

  const data = await response.json();

  // Store role in localStorage
  localStorage.setItem(USER_ROLE_KEY, data.role);

  return data;
}

/**
 * Get stored user role.
 */
export function getUserRole(): UserRole | null {
  return localStorage.getItem(USER_ROLE_KEY) as UserRole | null;
}

/**
 * Register current guest user with chosen username and password.
 *
 * @param userName - Desired username (unique, immutable)
 * @param password - Password to set
 * @returns New token and username
 */
export async function registerUser(userName: string, password: string): Promise<RegisterResponse> {
  const token = getToken();
  if (!token) {
    throw new Error('Not authenticated');
  }

  const response = await fetch('/api/user/register', {
    method: 'POST',
    headers: {
      'Authorization': `Bearer ${token}`,
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({ userName, password }),
  });

  if (!response.ok) {
    const errorData = await response.json().catch(() => ({ error: response.statusText }));
    throw new Error(errorData.error || `Registration failed: ${response.statusText}`);
  }

  const data = await response.json();

  // Update stored data with new token and username
  localStorage.setItem(TOKEN_KEY, data.token);
  localStorage.setItem(USERNAME_KEY, data.userName);
  localStorage.setItem(USER_ROLE_KEY, data.role);

  return data;
}

/**
 * Change current user's password.
 */
export async function changePassword(oldPassword: string, newPassword: string): Promise<void> {
  const token = getToken();
  if (!token) {
    throw new Error('Not authenticated');
  }

  const response = await fetch('/api/user/change-password', {
    method: 'POST',
    headers: {
      'Authorization': `Bearer ${token}`,
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({ oldPassword, newPassword }),
  });

  if (!response.ok) {
    const errorData = await response.json().catch(() => ({ error: response.statusText }));
    throw new Error(errorData.error || `Password change failed: ${response.statusText}`);
  }
}
