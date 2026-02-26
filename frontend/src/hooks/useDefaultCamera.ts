import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import {
	type DefaultCameraResponse,
	getDefaultCamera,
	setDefaultCamera,
} from "../myapi/client";
import { useAppStore } from "../store";

const defaultCameraKeys = {
	all: (roomId: string) => ["defaultCamera", roomId] as const,
};

export const useDefaultCamera = () => {
	const roomId = useAppStore((state) => state.roomId);
	const queryClient = useQueryClient();

	const query = useQuery({
		queryKey: defaultCameraKeys.all(roomId!),
		queryFn: () => getDefaultCamera(roomId!),
		enabled: !!roomId,
		staleTime: 30_000,
	});

	const mutation = useMutation({
		mutationFn: (cameraKey: string | null) =>
			setDefaultCamera(roomId!, cameraKey),
		onSuccess: (data) => {
			queryClient.setQueryData<DefaultCameraResponse>(
				defaultCameraKeys.all(roomId!),
				data,
			);
		},
	});

	return {
		defaultCamera: query.data?.default_camera ?? null,
		isLoading: query.isLoading,
		setDefaultCamera: mutation.mutate,
		isSettingDefault: mutation.isPending,
	};
};
