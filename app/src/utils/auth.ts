/**
 * JWT Authentication utilities for browser client.
 *
 * Uses jose library for secure JWT decoding and validation.
 * JWT is the single source of truth for user identity and role.
 *
 * Authentication Flow:
 * 1. Call /api/user/register to create user (guest or with password)
 * 2. Call /api/login to authenticate and get JWT token
 *
 * Architecture:
 * - userName is the primary identifier (unique, immutable after registration)
 * - Guests get auto-generated usernames (user-xyz)
 * - Upgrade allows choosing permanent username + password
 * - Role is always read from JWT (single source of truth)
 */

import { decodeJwt, type JWTPayload } from "jose";

const TOKEN_KEY = "zndraw_jwt_token";

export type UserRole = "guest" | "user" | "admin";

/** JWT payload structure from server */
interface ZnDrawJWTPayload extends JWTPayload {
	sub: string; // userName
	role: UserRole;
	jti: string; // JWT ID
	iat: number; // Issued at
	exp: number; // Expiration (7 days)
}

export interface LoginResponse {
	status: string;
	token: string;
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
 * For guests: First registers user via /api/user/register, then authenticates via /api/login.
 * For registered users with password: Authenticates directly via /api/login.
 *
 * @param userName - Optional username (if not provided, server generates guest username)
 * @param password - Optional password (if provided, authenticates as registered user)
 */
export async function login(
	userName?: string,
	password?: string,
): Promise<LoginResponse> {
	// For guests (no password), we need to register first, then login
	if (!password) {
		// Step 1: Register user (creates in backend)
		const registerBody: { userName?: string } = {};
		if (userName !== undefined) {
			registerBody.userName = userName;
		}

		const registerResponse = await fetch("/api/user/register", {
			method: "POST",
			headers: { "Content-Type": "application/json" },
			body: JSON.stringify(registerBody),
		});

		// Handle different registration outcomes
		let usernameToLogin: string;

		if (registerResponse.ok) {
			// New user created - get username from response
			const registerData = await registerResponse.json();
			usernameToLogin = registerData.userName;
		} else if (registerResponse.status === 409) {
			// User already exists - use provided username
			if (!userName) {
				// Should never happen: 409 without providing username
				throw new Error("Registration conflict without username");
			}
			usernameToLogin = userName;
		} else {
			// Other error - throw
			const errorData = await registerResponse
				.json()
				.catch(() => ({ error: registerResponse.statusText }));
			throw new Error(
				errorData.error ||
					`Registration failed: ${registerResponse.statusText}`,
			);
		}

		// Step 2: Login (get JWT)
		const loginResponse = await fetch("/api/login", {
			method: "POST",
			headers: { "Content-Type": "application/json" },
			body: JSON.stringify({ userName: usernameToLogin }),
		});

		if (!loginResponse.ok) {
			const errorData = await loginResponse
				.json()
				.catch(() => ({ error: loginResponse.statusText }));
			throw new Error(
				errorData.error || `Login failed: ${loginResponse.statusText}`,
			);
		}

		const data = await loginResponse.json();
		localStorage.setItem(TOKEN_KEY, data.token);
		return data;
	}

	// For registered users with password: direct login
	const response = await fetch("/api/login", {
		method: "POST",
		headers: { "Content-Type": "application/json" },
		body: JSON.stringify({ userName, password }),
	});

	if (!response.ok) {
		const errorData = await response
			.json()
			.catch(() => ({ error: response.statusText }));
		throw new Error(errorData.error || `Login failed: ${response.statusText}`);
	}

	const data = await response.json();
	localStorage.setItem(TOKEN_KEY, data.token);
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
 * Clear stored authentication data.
 */
export function logout(): void {
	localStorage.removeItem(TOKEN_KEY);
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
 * Check if stored token is valid (exists, can be decoded, and not expired).
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

	// Check expiry (tokens now expire after 7 days)
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
 * Get username from token claims.
 */
export function getUsernameFromToken(): string | null {
	const claims = getTokenClaims();
	return claims?.sub ?? null;
}

/**
 * Get user role from JWT token (single source of truth).
 */
export function getUserRole(): UserRole | null {
	const claims = getTokenClaims();
	return claims?.role ?? null;
}

/**
 * Upgrade current guest user to registered with chosen username and password.
 *
 * Requires authentication (must be logged in as guest first).
 * After upgrade, guest becomes a registered user with permanent credentials.
 *
 * @param userName - Desired username (unique, immutable)
 * @param password - Password to set (must meet minimum requirements)
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

	const response = await fetch("/api/user/upgrade", {
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
			errorData.error || `Upgrade failed: ${response.statusText}`,
		);
	}

	const data = await response.json();

	// Update stored token - role and username come from new JWT
	localStorage.setItem(TOKEN_KEY, data.token);

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
