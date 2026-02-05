import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import { getDefaultCamera, setDefaultCamera } from "../myapi/client";
import { useAppStore } from "../store";

/**
 * Hook for managing the room's default camera.
 *
 * The default camera is cloned to new session cameras when users join the room.
 */
export const useDefaultCamera = () => {
	const roomId = useAppStore((state) => state.roomId);
	const queryClient = useQueryClient();
	const showSnackbar = useAppStore((state) => state.showSnackbar);

	const { data, isLoading } = useQuery({
		queryKey: ["defaultCamera", roomId],
		queryFn: () => getDefaultCamera(roomId!),
		enabled: !!roomId,
		staleTime: 30000, // 30 seconds
	});

	const { mutate: setDefault, isPending: isSettingDefault } = useMutation({
		mutationFn: (cameraKey: string | null) =>
			setDefaultCamera(roomId!, cameraKey),
		onSuccess: (_, cameraKey) => {
			queryClient.invalidateQueries({ queryKey: ["defaultCamera", roomId] });
			if (cameraKey) {
				showSnackbar(
					`Set "${cameraKey}" as default camera for new sessions`,
					"success",
				);
			} else {
				showSnackbar("Cleared default camera", "info");
			}
		},
		onError: (error: any) => {
			const message =
				error?.response?.data?.error || "Failed to update default camera";
			showSnackbar(message, "error");
		},
	});

	return {
		defaultCamera: data?.default_camera ?? null,
		isLoading,
		isSettingDefault,
		setDefaultCamera: setDefault,
	};
};
