import SearchIcon from "@mui/icons-material/Search";
import {
	Alert,
	Box,
	Button,
	CircularProgress,
	Divider,
	InputAdornment,
	TextField,
} from "@mui/material";
import { useQuery } from "@tanstack/react-query";
import { isAxiosError } from "axios";
import React, { useState } from "react";
import { useNavigate } from "react-router-dom";
import {
	FileBreadcrumbs,
	FileList,
	FilesystemSelector,
	LoadFileDialog,
} from "../components/filesystem";
import type { LoadFileParams } from "../components/filesystem/LoadFileDialog";
import { useFilesystemProviders } from "../hooks/useFilesystemProviders";
import { useLeaveRoom } from "../hooks/useLeaveRoom";
import {
	createRoom,
	type FilesystemFileItem,
	listRooms,
	type ProviderInfo,
	readProvider,
	submitTask,
} from "../myapi/client";
import { useAppStore } from "../store";
import { useDockviewApi } from "../stores/dockviewApiStore";

/**
 * FilesystemPanel — sidebar-embedded filesystem browser.
 *
 * Renders inside the dockview sidebar zone (no AppBar / Container /
 * Paper chrome). Business logic (queries, state, handlers, dialog)
 * lives here directly.
 *
 * Notes
 * -----
 * - ``roomId`` comes from ``useAppStore`` so the panel works regardless
 *   of route.
 * - Post-load navigation runs the ``useLeaveRoom`` cascade before
 *   ``navigate`` so plots / viewer state don't leak across rooms.
 */
export function FilesystemPanel() {
	const navigate = useNavigate();
	const roomId = useAppStore((s) => s.roomId);
	const showSnackbar = useAppStore((state) => state.showSnackbar);
	const dockApi = useDockviewApi((s) => s.api);
	const leaveRoom = useLeaveRoom({ api: dockApi });

	const [selectedProvider, setSelectedProvider] = useState<ProviderInfo | null>(
		null,
	);
	const [currentPath, setCurrentPath] = useState<string>("/");
	const [globPattern, setGlobPattern] = useState<string>("");
	const [activeGlob, setActiveGlob] = useState<string>("");
	const [loadDialog, setLoadDialog] = useState<{
		open: boolean;
		file: FilesystemFileItem | null;
	}>({ open: false, file: null });
	const [isLoadingFile, setIsLoadingFile] = useState(false);

	const clearSearch = () => {
		setGlobPattern("");
		setActiveGlob("");
	};

	// Query for available filesystem providers
	const {
		data: providers,
		isLoading: isLoadingProviders,
		error: providersError,
	} = useFilesystemProviders();

	// Query for rooms list (used by "Existing room" dropdown in LoadFileDialog)
	const { data: rooms } = useQuery({
		queryKey: ["rooms"],
		queryFn: () => listRooms(),
	});

	// Auto-select first provider if only one is available
	React.useEffect(() => {
		if (providers && providers.length === 1 && !selectedProvider) {
			setSelectedProvider(providers[0]);
		}
	}, [providers, selectedProvider]);

	// Query for directory listing or glob search via provider read
	const providerParams: Record<string, string> = { path: currentPath };
	if (activeGlob) {
		providerParams.glob = activeGlob;
	}

	const {
		data: files,
		isLoading: isLoadingFiles,
		error: filesError,
	} = useQuery({
		queryKey: [
			"providerRead",
			roomId,
			selectedProvider?.full_name,
			currentPath,
			activeGlob,
		],
		queryFn: () => {
			if (!roomId || !selectedProvider) return Promise.resolve([]);
			return readProvider(roomId, selectedProvider.full_name, providerParams);
		},
		enabled: !!roomId && !!selectedProvider,
	});

	// Submit LoadFile task: resolve target room, cascade-close the current
	// room, navigate, then fire the async task in the background.
	const handleLoadFile = async (params: LoadFileParams) => {
		if (!selectedProvider || !roomId) return;

		setIsLoadingFile(true);

		const { room_target, ...taskParams } = params;

		// Resolve target room ID (create if needed)
		let targetRoomId: string;
		try {
			if (room_target.type === "new") {
				targetRoomId = room_target.room_id;
				await createRoom({
					room_id: targetRoomId,
					copy_from: "@none",
					...(room_target.description && {
						description: room_target.description,
					}),
				});
			} else if (room_target.type === "existing") {
				targetRoomId = room_target.room_id;
			} else {
				targetRoomId = roomId;
			}
		} catch (error: unknown) {
			const detail = isAxiosError(error)
				? (error.response?.data?.detail ?? error.message)
				: "Failed to create room";
			showSnackbar(detail, "error");
			setIsLoadingFile(false);
			return;
		}

		// Close the load dialog and cascade-close the current room so plots /
		// viewer state don't leak across rooms, then navigate to the target.
		setLoadDialog({ open: false, file: null });
		setIsLoadingFile(false);
		if (targetRoomId !== roomId) {
			await leaveRoom({ skipConfirm: true });
		}
		navigate(`/rooms/${targetRoomId}`);

		// Submit the task in the background (returns 202 ACCEPTED)
		const loadFileJobName = `${selectedProvider.room_id}:modifiers:LoadFile`;
		submitTask(targetRoomId, loadFileJobName, {
			provider_name: selectedProvider.full_name,
			...taskParams,
		}).then(
			() => showSnackbar("File loading task submitted", "success"),
			(error: unknown) => {
				const detail = isAxiosError(error)
					? (error.response?.data?.detail ?? error.message)
					: "Failed to load file";
				showSnackbar(detail, "error");
			},
		);
	};

	// Navigation handlers
	const handleItemClick = (item: FilesystemFileItem) => {
		if (item.type === "directory") {
			setCurrentPath(item.path);
			clearSearch();
		} else {
			setLoadDialog({ open: true, file: item });
		}
	};

	const handleNavigate = (path: string) => {
		setCurrentPath(path);
		clearSearch();
	};

	const handleGoToRoot = () => {
		setCurrentPath("/");
		clearSearch();
	};

	const handleGoToParent = () => {
		const pathParts = currentPath.split("/").filter(Boolean);
		if (pathParts.length > 1) {
			setCurrentPath(`/${pathParts.slice(0, -1).join("/")}`);
		} else {
			setCurrentPath("/");
		}
		clearSearch();
	};

	const handleProviderSelect = (fullName: string) => {
		const provider = providers?.find((p) => p.full_name === fullName);
		if (provider) {
			setSelectedProvider(provider);
			setCurrentPath("/");
			clearSearch();
		}
	};

	const handleSearchSubmit = (e: React.FormEvent) => {
		e.preventDefault();
		setActiveGlob(globPattern);
	};

	// Show error state for providers
	if (providersError) {
		return (
			<Box
				data-testid="filesystem-panel"
				sx={{
					display: "flex",
					flexDirection: "column",
					height: "100%",
					overflow: "auto",
					p: 2,
				}}
			>
				<Alert severity="error">
					{isAxiosError(providersError)
						? (providersError.response?.data?.detail ?? providersError.message)
						: "Failed to load filesystem providers"}
				</Alert>
			</Box>
		);
	}

	// Show loading state
	if (isLoadingProviders) {
		return (
			<Box
				data-testid="filesystem-panel"
				sx={{
					display: "flex",
					flexDirection: "column",
					height: "100%",
					overflow: "auto",
					p: 2,
				}}
			>
				<Box sx={{ display: "flex", justifyContent: "center", p: 4 }}>
					<CircularProgress />
				</Box>
			</Box>
		);
	}

	// Show no providers message
	if (providers && providers.length === 0) {
		return (
			<Box
				data-testid="filesystem-panel"
				sx={{
					display: "flex",
					flexDirection: "column",
					height: "100%",
					overflow: "auto",
					p: 2,
				}}
			>
				<Alert severity="info">
					No filesystems are currently registered for this room. Register a
					filesystem using the ZnDraw Python client.
				</Alert>
			</Box>
		);
	}

	return (
		<Box
			data-testid="filesystem-panel"
			sx={{
				display: "flex",
				flexDirection: "column",
				height: "100%",
				overflow: "auto",
				p: 2,
			}}
		>
			{/* Provider Selector */}
			{providers && (
				<FilesystemSelector
					providers={providers}
					selected={selectedProvider?.full_name ?? ""}
					onSelect={handleProviderSelect}
				/>
			)}

			{/* Breadcrumbs + Search */}
			{selectedProvider && (
				<>
					<FileBreadcrumbs
						currentPath={currentPath}
						onNavigate={handleNavigate}
						onGoToRoot={handleGoToRoot}
					/>
					<Box component="form" onSubmit={handleSearchSubmit} sx={{ mb: 2 }}>
						<TextField
							size="small"
							fullWidth
							placeholder="Glob pattern (e.g. **/*.xyz)"
							value={globPattern}
							onChange={(e) => setGlobPattern(e.target.value)}
							slotProps={{
								input: {
									startAdornment: (
										<InputAdornment position="start">
											<SearchIcon />
										</InputAdornment>
									),
									endAdornment: activeGlob ? (
										<InputAdornment position="end">
											<Button size="small" onClick={clearSearch}>
												Clear
											</Button>
										</InputAdornment>
									) : null,
								},
							}}
						/>
					</Box>
					<Divider sx={{ mb: 2 }} />
				</>
			)}

			{/* File List */}
			<FileList
				files={files ?? []}
				currentPath={currentPath}
				isLoading={isLoadingFiles}
				error={filesError as Error | null}
				onItemClick={handleItemClick}
				onGoToParent={handleGoToParent}
			/>

			{/* Load File Dialog */}
			<LoadFileDialog
				open={loadDialog.open}
				file={loadDialog.file}
				isLoading={isLoadingFile}
				rooms={rooms ?? []}
				onClose={() => setLoadDialog({ open: false, file: null })}
				onLoad={handleLoadFile}
			/>
		</Box>
	);
}
