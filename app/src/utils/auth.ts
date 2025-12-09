/**
 * JWT Authentication utilities for browser client.
 *
 * Uses jose library for secure JWT decoding and validation.
 *
 * Simplified Architecture:
 * - userName is the primary identifier (unique, immutable after registration)
 * - Guests get auto-generated usernames (user-xyz)
 * - Registration allows choosing permanent username
 * - No client_id concept
 */

import { decodeJwt, type JWTPayload } from "jose";

const TOKEN_KEY = "zndraw_jwt_token";
const USERNAME_KEY = "zndraw_username";
const USER_ROLE_KEY = "zndraw_user_role";

export type UserRole = "guest" | "user" | "admin";

/** JWT payload structure from server */
interface ZnDrawJWTPayload extends JWTPayload {
	sub: string; // userName
	role: UserRole;
	jti: string; // JWT ID
}

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
export async function login(
	userName?: string,
	password?: string,
): Promise<LoginResponse> {
	const body: { userName?: string; password?: string } = {};

	if (userName !== undefined) {
		body.userName = userName;
	}
	if (password !== undefined) {
		body.password = password;
	}

	const response = await fetch("/api/login", {
		method: "POST",
		headers: { "Content-Type": "application/json" },
		body: JSON.stringify(body),
	});

	if (!response.ok) {
		const errorData = await response
			.json()
			.catch(() => ({ error: response.statusText }));
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
 * Ensure user is authenticated with a valid token. Auto-login if needed.
 * Call this before making authenticated requests.
 */
export async function ensureAuthenticated(): Promise<void> {
	if (!isTokenValid()) {
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
 * Check if user is authenticated (has valid token).
 */
export function isAuthenticated(): boolean {
	return isTokenValid();
}

/**
 * Decode JWT token and extract payload using jose library.
 *
 * @param token - JWT token string
 * @returns Decoded payload or null if invalid
 */
export function decodeToken(token: string): ZnDrawJWTPayload | null {
	try {
		return decodeJwt(token) as ZnDrawJWTPayload;
	} catch {
		return null;
	}
}

/**
 * Check if stored token is valid (exists and can be decoded).
 * If token has expiry, also checks if not expired.
 */
export function isTokenValid(): boolean {
	const token = getToken();
	if (!token) {
		return false;
	}

	const payload = decodeToken(token);
	if (!payload) {
		logout(); // Clear invalid token
		return false;
	}

	// Check expiry if present (currently tokens don't expire, but future-proof)
	if (payload.exp && Date.now() >= payload.exp * 1000) {
		logout(); // Clear expired token
		return false;
	}

	return true;
}

/**
 * Get claims from stored token without server request.
 *
 * @returns Token claims or null if no valid token
 */
export function getTokenClaims(): ZnDrawJWTPayload | null {
	const token = getToken();
	if (!token) {
		return null;
	}
	return decodeToken(token);
}

/**
 * Get username from token claims (preferred over localStorage).
 */
export function getUsernameFromToken(): string | null {
	const claims = getTokenClaims();
	return claims?.sub ?? null;
}

/**
 * Get role from token claims (preferred over localStorage).
 */
export function getRoleFromToken(): UserRole | null {
	const claims = getTokenClaims();
	return claims?.role ?? null;
}

/**
 * Fetch user role from server.
 */
export async function fetchUserRole(): Promise<UserRoleResponse> {
	const token = getToken();
	if (!token) {
		throw new Error("Not authenticated");
	}

	const response = await fetch("/api/user/role", {
		method: "GET",
		headers: {
			Authorization: `Bearer ${token}`,
			"Content-Type": "application/json",
		},
	});

	if (!response.ok) {
		const errorData = await response
			.json()
			.catch(() => ({ error: response.statusText }));
		throw new Error(
			errorData.error || `Failed to fetch role: ${response.statusText}`,
		);
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
export async function registerUser(
	userName: string,
	password: string,
): Promise<RegisterResponse> {
	const token = getToken();
	if (!token) {
		throw new Error("Not authenticated");
	}

	const response = await fetch("/api/user/register", {
		method: "POST",
		headers: {
			Authorization: `Bearer ${token}`,
			"Content-Type": "application/json",
		},
		body: JSON.stringify({ userName, password }),
	});

	if (!response.ok) {
		const errorData = await response
			.json()
			.catch(() => ({ error: response.statusText }));
		throw new Error(
			errorData.error || `Registration failed: ${response.statusText}`,
		);
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
export async function changePassword(
	oldPassword: string,
	newPassword: string,
): Promise<void> {
	const token = getToken();
	if (!token) {
		throw new Error("Not authenticated");
	}

	const response = await fetch("/api/user/change-password", {
		method: "POST",
		headers: {
			Authorization: `Bearer ${token}`,
			"Content-Type": "application/json",
		},
		body: JSON.stringify({ oldPassword, newPassword }),
	});

	if (!response.ok) {
		const errorData = await response
			.json()
			.catch(() => ({ error: response.statusText }));
		throw new Error(
			errorData.error || `Password change failed: ${response.statusText}`,
		);
	}
}
