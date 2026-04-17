import { useQuery, useQueryClient } from "@tanstack/react-query";
import { useEffect } from "react";
import { listProviders } from "../myapi/client";
import { socket } from "../socket";
import { useAppStore } from "../store";

/**
 * Returns ``true`` when the current room has at least one provider of
 * category ``"filesystem"``. Invalidates on the ``providers_invalidate``
 * Socket.IO event so a newly-registered provider lights the icon up
 * without a page reload.
 */
export function useHasFilesystemProviders(): boolean {
	const roomId = useAppStore((s) => s.roomId);
	const queryClient = useQueryClient();

	const { data } = useQuery({
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

	return (data?.length ?? 0) > 0;
}
