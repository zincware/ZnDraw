import ArrowBackIcon from "@mui/icons-material/ArrowBack";
import FolderOpenIcon from "@mui/icons-material/FolderOpen";
import SearchIcon from "@mui/icons-material/Search";
import {
	Alert,
	AppBar,
	Box,
	Button,
	CircularProgress,
	Container,
	Divider,
	IconButton,
	InputAdornment,
	Paper,
	TextField,
	Toolbar,
	Typography,
} from "@mui/material";
import { useQuery } from "@tanstack/react-query";
import { isAxiosError } from "axios";
import React, { useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import {
	FileBreadcrumbs,
	FileList,
	FilesystemSelector,
	LoadFileDialog,
} from "../components/filesystem";
import type { LoadFileParams } from "../components/filesystem/LoadFileDialog";
import {
	type FilesystemFileItem,
	type ProviderInfo,
	createRoom,
	listProviders,
	listRooms,
	readProvider,
	submitTask,
} from "../myapi/client";
import { useAppStore } from "../store";

/**
 * Filesystem browser page for browsing registered filesystem providers
 * and loading files into ZnDraw rooms via the joblib task system.
 */
export default function FilesystemBrowserPage() {
	const navigate = useNavigate();
	const { roomId } = useParams<{ roomId: string }>();
	const showSnackbar = useAppStore((state) => state.showSnackbar);

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
	} = useQuery({
		queryKey: ["filesystemProviders", roomId],
		queryFn: () => listProviders(roomId!, "filesystem"),
		enabled: !!roomId,
		retry: false,
	});

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
		queryFn: () =>
			readProvider(roomId!, selectedProvider!.full_name, providerParams),
		enabled: !!roomId && !!selectedProvider,
	});

	// Submit LoadFile task: resolve target room, navigate immediately,
	// then fire the async task in the background.
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

		// Navigate to the target room immediately
		setLoadDialog({ open: false, file: null });
		setIsLoadingFile(false);
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
			setCurrentPath("/" + pathParts.slice(0, -1).join("/"));
		} else {
			setCurrentPath("/");
		}
		clearSearch();
	};

	const handleBack = () => {
		if (roomId) {
			navigate(`/rooms/${roomId}`);
		} else {
			navigate("/");
		}
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
			<Container maxWidth="lg" sx={{ mt: 4 }}>
				<Alert severity="error">
					{isAxiosError(providersError) ? (providersError.response?.data?.detail ?? providersError.message) : "Failed to load filesystem providers"}
				</Alert>
				<Button onClick={handleBack} sx={{ mt: 2 }}>
					Back
				</Button>
			</Container>
		);
	}

	// Show loading state
	if (isLoadingProviders) {
		return (
			<Container maxWidth="lg" sx={{ mt: 4 }}>
				<Box sx={{ display: "flex", justifyContent: "center", p: 4 }}>
					<CircularProgress />
				</Box>
			</Container>
		);
	}

	// Show no providers message
	if (providers && providers.length === 0) {
		return (
			<Container maxWidth="lg" sx={{ mt: 4 }}>
				<Alert severity="info">
					No filesystems are currently registered for this room. Register a
					filesystem using the ZnDraw Python client.
				</Alert>
				<Button onClick={handleBack} sx={{ mt: 2 }}>
					Back
				</Button>
			</Container>
		);
	}

	return (
		<Box sx={{ flexGrow: 1 }}>
			<AppBar position="static">
				<Toolbar>
					<IconButton
						edge="start"
						color="inherit"
						onClick={handleBack}
						sx={{ mr: 2 }}
					>
						<ArrowBackIcon />
					</IconButton>
					<FolderOpenIcon sx={{ mr: 1 }} />
					<Typography variant="h6" component="div" sx={{ flexGrow: 1 }}>
						Filesystem
					</Typography>
				</Toolbar>
			</AppBar>

			<Container maxWidth="lg" sx={{ mt: 4, mb: 4 }}>
				<Paper sx={{ p: 3 }}>
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
							<Box
								component="form"
								onSubmit={handleSearchSubmit}
								sx={{ mb: 2 }}
							>
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
				</Paper>
			</Container>

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
