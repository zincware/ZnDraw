import { useEffect, useCallback, useRef } from "react";
import { useAppStore } from "../store";
import { useParams } from "react-router-dom";
import {
	ensureAuthenticated,
	getUserRole,
	getUsernameFromToken,
} from "../utils/auth";

/**
 * Hook to handle REST-based initialization before socket connection.
 *
 * This hook only handles authentication and sets the roomId.
 * Room join and data fetching is handled by useSocketManager via room:join.
 */
export const useRestJoinManager = () => {
	const { setRoomId, setUserName, setUserRole } = useAppStore();
	const { roomId: room } = useParams<{ roomId: string }>();
	const hasInitialized = useRef(false);

	const initialize = useCallback(async () => {
		if (!room || hasInitialized.current) {
			return;
		}

		try {
			// Ensure user is authenticated before socket connection
			await ensureAuthenticated();

			// Get user role from JWT (single source of truth)
			const userRole = getUserRole();
			if (userRole) {
				setUserRole(userRole);
			}

			// Get username from auth system
			const userName = getUsernameFromToken();
			if (userName) {
				setUserName(userName);
			}

			// Set roomId to trigger socket connection
			setRoomId(room);
			hasInitialized.current = true;
		} catch (error) {
			console.error("Authentication failed:", error);
		}
	}, [room, setRoomId, setUserName, setUserRole]);

	useEffect(() => {
		initialize();

		return () => {
			// Reset on unmount for clean re-initialization
			hasInitialized.current = false;
		};
	}, [initialize]);
};
