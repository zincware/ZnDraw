import { useQuery, useQueryClient } from "@tanstack/react-query";
import { useEffect } from "react";
import { listProviders } from "../myapi/client";
import { socket } from "../socket";
import { useAppStore } from "../store";

/**
 * Query hook for filesystem providers in the current room.
 *
 * Sets ``staleTime: 5_000`` and subscribes to the
 * ``providers_invalidate`` Socket.IO event so a newly-registered
 * provider triggers a refetch without a page reload. Shared by
 * ``useHasFilesystemProviders`` and ``FilesystemPanel`` so both
 * consumers agree on cache policy and invalidation.
 */
export function useFilesystemProviders() {
	const roomId = useAppStore((s) => s.roomId);
	const queryClient = useQueryClient();

	const result = useQuery({
		queryKey: ["filesystemProviders", roomId],
		queryFn: () => listProviders(roomId!, "filesystem"),
		enabled: !!roomId,
		retry: false,
		staleTime: 5_000,
	});

	useEffect(() => {
		if (!roomId) return;
		const handle = () => {
			queryClient.invalidateQueries({
				queryKey: ["filesystemProviders", roomId],
			});
		};
		socket.on("providers_invalidate", handle);
		return () => {
			socket.off("providers_invalidate", handle);
		};
	}, [roomId, queryClient]);

	return result;
}
