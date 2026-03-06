/**
 * Authentication utilities for browser client.
 *
 * Uses fastapi-users backend (zndraw-auth):
 * - Guest: POST /v1/auth/guest → {access_token, token_type}
 * - Login: POST /v1/auth/jwt/login → {access_token, token_type}
 * - Register: POST /v1/auth/register → UserRead (no token)
 * - User info: GET /v1/auth/users/me → UserRead
 * - Change password: PATCH /v1/auth/users/me → UserRead
 */

const TOKEN_KEY = "zndraw_jwt_token";

export type UserRole = "user" | "admin";

export interface UserInfo {
	id: string;
	email: string;
	is_active: boolean;
	is_superuser: boolean;
	is_verified: boolean;
}

export interface AuthResult {
	token: string;
	user: UserInfo;
}

/**
 * Fetch current user info using a specific token.
 */
async function fetchMeWithToken(token: string | null): Promise<UserInfo> {
	if (!token) {
		throw new Error("No token provided");
	}
	const response = await fetch("/v1/auth/users/me", {
		headers: { Authorization: `Bearer ${token}` },
	});
	if (!response.ok) {
		throw new Error(`Failed to fetch user info: ${response.statusText}`);
	}
	return response.json();
}

/**
 * Login and get JWT token from server.
 *
 * For guests: POST /v1/auth/guest.
 * For registered users: POST /v1/auth/jwt/login with OAuth2 form data.
 *
 * @param email - Optional email (required for password login)
 * @param password - Optional password (if provided, authenticates as registered user)
 */
export async function login(
	email?: string,
	password?: string,
): Promise<AuthResult> {
	if (!password) {
		// Guest login
		const response = await fetch("/v1/auth/guest", { method: "POST" });
		if (!response.ok) {
			const errorData = await response
				.json()
				.catch(() => ({ detail: response.statusText }));
			throw new Error(
				errorData.detail || `Guest login failed: ${response.statusText}`,
			);
		}
		const data = await response.json();
		// Validate token before storing — prevents stale token on fetchMe failure
		const user = await fetchMeWithToken(data.access_token);
		localStorage.setItem(TOKEN_KEY, data.access_token);
		return { token: data.access_token, user };
	}

	// Password login (OAuth2 form data)
	const formData = new URLSearchParams();
	formData.append("username", email || "");
	formData.append("password", password);

	const response = await fetch("/v1/auth/jwt/login", {
		method: "POST",
		headers: { "Content-Type": "application/x-www-form-urlencoded" },
		body: formData,
	});
	if (!response.ok) {
		const errorData = await response
			.json()
			.catch(() => ({ detail: response.statusText }));
		throw new Error(errorData.detail || `Login failed: ${response.statusText}`);
	}
	const data = await response.json();
	// Validate token before storing — prevents stale token on fetchMe failure
	const user = await fetchMeWithToken(data.access_token);
	localStorage.setItem(TOKEN_KEY, data.access_token);
	return { token: data.access_token, user };
}

/**
 * Register a new user with email and password.
 *
 * POST /v1/auth/register returns UserRead (no token), so auto-login after.
 */
export async function registerUser(
	email: string,
	password: string,
): Promise<AuthResult> {
	const response = await fetch("/v1/auth/register", {
		method: "POST",
		headers: { "Content-Type": "application/json" },
		body: JSON.stringify({ email, password }),
	});
	if (!response.ok) {
		const errorData = await response
			.json()
			.catch(() => ({ detail: response.statusText }));
		throw new Error(
			errorData.detail || `Registration failed: ${response.statusText}`,
		);
	}
	// Register returns UserRead but no token — auto-login to get JWT
	return login(email, password);
}

/**
 * Change current user's password via PATCH /v1/auth/users/me.
 */
export async function changePassword(newPassword: string): Promise<void> {
	const token = getToken();
	const response = await fetch("/v1/auth/users/me", {
		method: "PATCH",
		headers: {
			"Content-Type": "application/json",
			Authorization: `Bearer ${token}`,
		},
		body: JSON.stringify({ password: newPassword }),
	});
	if (!response.ok) {
		const errorData = await response
			.json()
			.catch(() => ({ detail: response.statusText }));
		throw new Error(
			errorData.detail || `Password change failed: ${response.statusText}`,
		);
	}
}

/**
 * Ensure a valid token exists in localStorage, obtaining a fresh guest token if needed.
 *
 * Validates existing tokens server-side before returning. Coalesces concurrent
 * calls to prevent duplicate guest user creation (e.g., 401 interceptor + socket
 * flow racing on page load).
 */
let acquirePromise: Promise<AuthResult> | null = null;

export function acquireToken(): Promise<AuthResult> {
	if (acquirePromise) return acquirePromise;
	acquirePromise = doAcquireToken().finally(() => {
		acquirePromise = null;
	});
	return acquirePromise;
}

async function doAcquireToken(): Promise<AuthResult> {
	const existing = getToken();
	if (existing) {
		try {
			const user = await fetchMeWithToken(existing);
			return { token: existing, user };
		} catch {
			logout();
		}
	}
	return login();
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
