import AccountCircleIcon from "@mui/icons-material/AccountCircle";
import AdminPanelSettingsIcon from "@mui/icons-material/AdminPanelSettings";
import BrushIcon from "@mui/icons-material/Brush";
import CameraAltIcon from "@mui/icons-material/CameraAlt";
import CodeIcon from "@mui/icons-material/Code";
import DarkModeIcon from "@mui/icons-material/DarkMode";
import LightModeIcon from "@mui/icons-material/LightMode";
import LoginIcon from "@mui/icons-material/Login";
import LogoutIcon from "@mui/icons-material/Logout";
import ManageAccountsIcon from "@mui/icons-material/ManageAccounts";
import OpenWithIcon from "@mui/icons-material/OpenWith";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import Alert from "@mui/material/Alert";
import AppBar from "@mui/material/AppBar";
import Box from "@mui/material/Box";
import CssBaseline from "@mui/material/CssBaseline";
import IconButton from "@mui/material/IconButton";
import Link from "@mui/material/Link";
import ListItemIcon from "@mui/material/ListItemIcon";
import ListItemText from "@mui/material/ListItemText";
import Menu from "@mui/material/Menu";
import MenuItem from "@mui/material/MenuItem";
import Snackbar from "@mui/material/Snackbar";
import { useColorScheme } from "@mui/material/styles";
import Toolbar from "@mui/material/Toolbar";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import { useQueryClient } from "@tanstack/react-query";
import React, { useCallback, useEffect, useState } from "react";
import { useParams, useSearchParams } from "react-router-dom";
import AdminPanel from "../components/AdminPanel";
import ConnectionDialog from "../components/ConnectionDialog";
import LoginDialog from "../components/LoginDialog";
import FrameProgressBar from "../components/ProgressBar";
import ProgressNotifications from "../components/ProgressNotifications";
import RegisterDialog from "../components/RegisterDialog";
import RoomManagementMenu from "../components/RoomManagementMenu";
import SiMGenButtons from "../components/SiMGenButtons";
import SiMGenTutorialDialog from "../components/SiMGenTutorialDialog";
import UserProfileDialog from "../components/UserProfileDialog";
import { TRAJECTORY_ACCEPT } from "../constants/fileTypes";
import { LAYOUT_CONSTANTS } from "../constants/layout";
import { useKeyboardShortcuts } from "../hooks/useKeyboardShortcuts";
import { useSocketManager } from "../hooks/useSocketManager";
import { uploadTrajectory } from "../myapi/client";
import {
	ActivityBar,
	BottomZone,
	DockviewLayout,
	ensureViewerPanel,
	PANELS,
	type PanelId,
	resetDockview,
	SidebarZone,
} from "../panels";
import { useDockviewApi } from "../stores/dockviewApiStore";
import { connectWithAuth } from "../socket";
import { useAppStore } from "../store";
import { logout as authLogout } from "../utils/auth";
import { downloadScreenshot } from "../utils/screenshot";

export default function MainPage() {
	const { roomId } = useParams<{ roomId: string }>();
	const setRoomId = useAppStore((state) => state.setRoomId);

	// Set roomId in store for child components that read from it
	useEffect(() => {
		if (roomId) setRoomId(roomId);
	}, [roomId, setRoomId]);

	useSocketManager({ roomId });
	useKeyboardShortcuts();

	// Single source of truth for panel-drag lifecycle. Uses window-level
	// dragstart/dragend/drop plus a document keydown fallback so stuck
	// "hot" zones can't survive an Escape-cancel or a missed dragend.
	useEffect(() => {
		const DRAG_MIME = "application/x-zndraw-panel-id";
		const setActive = useAppStore.getState().setPanelDragActive;
		const onStart = (e: DragEvent) => {
			if (e.dataTransfer?.types.includes(DRAG_MIME)) setActive(true);
		};
		const onEnd = () => setActive(false);
		const onKey = (e: KeyboardEvent) => {
			if (e.key === "Escape") setActive(false);
		};
		window.addEventListener("dragstart", onStart, { passive: true });
		window.addEventListener("dragend", onEnd, { passive: true });
		window.addEventListener("drop", onEnd, { passive: true });
		document.addEventListener("keydown", onKey);
		return () => {
			window.removeEventListener("dragstart", onStart);
			window.removeEventListener("dragend", onEnd);
			window.removeEventListener("drop", onEnd);
			document.removeEventListener("keydown", onKey);
		};
	}, []);

	const interactionMode = useAppStore((state) => state.mode);
	const enterDrawingMode = useAppStore((state) => state.enterDrawingMode);
	const exitDrawingMode = useAppStore((state) => state.exitDrawingMode);
	const enterEditingMode = useAppStore((state) => state.enterEditingMode);
	const exitEditingMode = useAppStore((state) => state.exitEditingMode);
	const serverVersion = useAppStore((state) => state.serverVersion);
	const globalSettings = useAppStore((state) => state.globalSettings);
	const userEmail = useAppStore((state) => state.user?.email ?? null);
	const isAdmin = useAppStore((state) => state.user?.is_superuser ?? false);
	const setUser = useAppStore((state) => state.setUser);
	const showSnackbar = useAppStore((state) => state.showSnackbar);
	const snackbar = useAppStore((state) => state.snackbar);
	const hideSnackbar = useAppStore((state) => state.hideSnackbar);
	const screenshotCapture = useAppStore((state) => state.screenshotCapture);
	const activeLeft = useAppStore((state) => state.activeLeft);
	const activeRight = useAppStore((state) => state.activeRight);
	const activeBottom = useAppStore((state) => state.activeBottom);
	const toggleActive = useAppStore((state) => state.toggleActive);
	const [connectionDialogOpen, setConnectionDialogOpen] = useState(false);
	const [loginDialogOpen, setLoginDialogOpen] = useState(false);
	const [registerDialogOpen, setRegisterDialogOpen] = useState(false);
	const [adminPanelOpen, setAdminPanelOpen] = useState(false);
	const [userProfileDialogOpen, setUserProfileDialogOpen] = useState(false);
	const [profileAnchorEl, setProfileAnchorEl] = useState<null | HTMLElement>(
		null,
	);
	const profileMenuOpen = Boolean(profileAnchorEl);
	const { mode: colorMode, setMode: setColorMode } = useColorScheme();
	const queryClient = useQueryClient();

	// Drag and drop file upload
	const handleFiles = useCallback(
		async (files: File[]) => {
			if (!roomId) {
				showSnackbar("Cannot upload: no room selected", "warning");
				return;
			}
			for (const file of files) {
				try {
					await uploadTrajectory(roomId, file);
					queryClient.invalidateQueries({ queryKey: ["frame", roomId] });
					showSnackbar(`Uploaded: ${file.name}`, "success");
				} catch (error: any) {
					const detail =
						error?.response?.data?.detail || `Upload failed: ${file.name}`;
					showSnackbar(detail, "error");
				}
			}
		},
		[roomId, queryClient, showSnackbar],
	);

	// Tutorial dialog state
	const [tutorialDialogOpen, setTutorialDialogOpen] = useState(false);

	// File upload ref for button click
	const fileInputRef = React.useRef<HTMLInputElement>(null);

	const handleTakeScreenshot = async () => {
		if (!screenshotCapture) {
			showSnackbar("Screenshot capture not available", "error");
			return;
		}

		try {
			const blob = await screenshotCapture();
			downloadScreenshot(blob);
			showSnackbar("Screenshot downloaded", "success");
		} catch (error) {
			showSnackbar(
				`Screenshot failed: ${error instanceof Error ? error.message : "Unknown error"}`,
				"error",
			);
		}
	};

	const handleFileUploadClick = () => {
		fileInputRef.current?.click();
	};

	const handleFileInputChange = async (
		event: React.ChangeEvent<HTMLInputElement>,
	) => {
		const files = event.target.files;
		if (!files || files.length === 0) return;

		await handleFiles(Array.from(files));

		// Reset input
		if (fileInputRef.current) {
			fileInputRef.current.value = "";
		}
	};

	const handleToggleColorMode = () => {
		setColorMode(colorMode === "light" ? "dark" : "light");
	};

	const handleProfileClick = (event: React.MouseEvent<HTMLElement>) => {
		setProfileAnchorEl(event.currentTarget);
	};

	const handleProfileClose = () => {
		setProfileAnchorEl(null);
	};

	const handleSwitchUser = () => {
		handleProfileClose();
		setLoginDialogOpen(true);
	};

	const handleLogout = async () => {
		handleProfileClose();

		try {
			authLogout();
			const { user } = await connectWithAuth();
			setUser(user);

			showSnackbar("Logged out, reconnected as guest", "info");
		} catch (err) {
			showSnackbar("Logout failed", "error");
		}
	};

	// Re-add the viewer panel on roomId change if it has been closed.
	useEffect(() => {
		const api = useDockviewApi.getState().api;
		if (!api) return;
		if (!roomId) return;
		ensureViewerPanel(api);
	}, [roomId]);

	// Auto-open a sidebar panel from `?panel=` query param (for redirect targets
	// like /rooms/:id/files → /rooms/:id?panel=filesystem).
	const [searchParams] = useSearchParams();
	useEffect(() => {
		const panel = searchParams.get("panel") as PanelId | null;
		if (!panel || !(panel in PANELS)) return;
		const def = PANELS[panel];
		if (def.kind !== "tool" || def.default.bar === "editor") return;
		const current =
			def.default.bar === "left"
				? activeLeft
				: def.default.bar === "right"
					? activeRight
					: activeBottom;
		if (current !== panel) {
			toggleActive(def.default.bar, panel);
		}
	}, [searchParams, activeLeft, activeRight, activeBottom, toggleActive]);

	return (
		<>
			<Box
				sx={{
					display: "flex",
					flexDirection: "column",
					height: "100vh",
					width: "100vw",
					overflow: "hidden",
				}}
			>
				<CssBaseline />

				{/* Header / AppBar */}
				<AppBar
					position="static"
					sx={{
						zIndex: 0,
						height: LAYOUT_CONSTANTS.APPBAR_HEIGHT,
					}}
				>
					<Toolbar
						sx={{
							minHeight: `${LAYOUT_CONSTANTS.APPBAR_HEIGHT}px !important`,
							height: LAYOUT_CONSTANTS.APPBAR_HEIGHT,
							padding: "0 16px",
						}}
					>
						<Tooltip title="View on GitHub">
							<Link
								href="https://github.com/zincware/ZnDraw"
								target="_blank"
								rel="noopener noreferrer"
								color="inherit"
								underline="none"
								sx={{ "&:hover": { textDecoration: "none" }, cursor: "pointer" }}
							>
								<Typography variant="h6" noWrap component="div">
									{serverVersion ? `ZnDraw ${serverVersion}` : "ZnDraw"}
									{globalSettings?.simgen?.enabled && (
										<>
											{" + "}
											<Link
												href="https://github.com/RokasEl/simgen"
												target="_blank"
												rel="noopener noreferrer"
												color="inherit"
												underline="hover"
												sx={{ cursor: "pointer" }}
											>
												SiMGen
											</Link>
										</>
									)}
								</Typography>
							</Link>
						</Tooltip>
						<Box sx={{ flexGrow: 1 }} />
						<SiMGenButtons
							onTutorialClick={() => setTutorialDialogOpen(true)}
						/>
						<Box sx={{ flexGrow: 1 }} />
						<Tooltip
							title={
								interactionMode === "drawing"
									? "Disable drawing mode (X)"
									: "Enable drawing mode (X)"
							}
						>
							<IconButton
								color="inherit"
								aria-label="toggle drawing mode"
								onClick={() =>
									interactionMode === "drawing"
										? exitDrawingMode()
										: enterDrawingMode(queryClient)
								}
								sx={{
									backgroundColor:
										interactionMode === "drawing"
											? "rgba(255, 255, 255, 0.2)"
											: "transparent",
								}}
							>
								<BrushIcon />
							</IconButton>
						</Tooltip>
						<Tooltip
							title={
								interactionMode === "editing"
									? "Exit editing mode (E)"
									: "Enter editing mode (E)"
							}
						>
							<IconButton
								color="inherit"
								aria-label="toggle editing mode"
								onClick={() =>
									interactionMode === "editing"
										? exitEditingMode()
										: enterEditingMode()
								}
								sx={{
									backgroundColor:
										interactionMode === "editing"
											? "rgba(255, 255, 255, 0.2)"
											: "transparent",
								}}
							>
								<OpenWithIcon />
							</IconButton>
						</Tooltip>
						<Tooltip
							title={
								colorMode === "light"
									? "Switch to dark mode"
									: "Switch to light mode"
							}
						>
							<IconButton
								color="inherit"
								aria-label="toggle theme"
								onClick={handleToggleColorMode}
							>
								{colorMode === "light" ? <DarkModeIcon /> : <LightModeIcon />}
							</IconButton>
						</Tooltip>
						<Tooltip title={"Python code connection info"}>
							<IconButton
								color="inherit"
								aria-label="show connection info"
								onClick={() => setConnectionDialogOpen(true)}
							>
								<CodeIcon />
							</IconButton>
						</Tooltip>
						<Tooltip title="Take screenshot">
							<IconButton
								color="inherit"
								aria-label="take screenshot"
								onClick={handleTakeScreenshot}
								disabled={!screenshotCapture}
							>
								<CameraAltIcon />
							</IconButton>
						</Tooltip>
						<Tooltip title="Upload file">
							<IconButton
								color="inherit"
								aria-label="upload file"
								onClick={handleFileUploadClick}
							>
								<UploadFileIcon />
							</IconButton>
						</Tooltip>
						<RoomManagementMenu />
						<Tooltip title="User profile">
							<IconButton
								color="inherit"
								aria-label="user profile"
								onClick={handleProfileClick}
							>
								<AccountCircleIcon />
							</IconButton>
						</Tooltip>
					</Toolbar>
				</AppBar>

				{/* Profile Menu */}
				<Menu
					anchorEl={profileAnchorEl}
					open={profileMenuOpen}
					onClose={handleProfileClose}
					anchorOrigin={{
						vertical: "bottom",
						horizontal: "right",
					}}
					transformOrigin={{
						vertical: "top",
						horizontal: "right",
					}}
				>
					<MenuItem disabled>
						<ListItemText
							primary={userEmail || "Unknown"}
							secondary={isAdmin ? "Admin" : "User"}
						/>
					</MenuItem>

					<MenuItem
						onClick={() => {
							handleProfileClose();
							setRegisterDialogOpen(true);
						}}
					>
						<ListItemIcon>
							<ManageAccountsIcon fontSize="small" />
						</ListItemIcon>
						<ListItemText>Register</ListItemText>
					</MenuItem>
					<MenuItem onClick={handleSwitchUser}>
						<ListItemIcon>
							<LoginIcon fontSize="small" />
						</ListItemIcon>
						<ListItemText>Login</ListItemText>
					</MenuItem>
					<MenuItem
						onClick={() => {
							handleProfileClose();
							setUserProfileDialogOpen(true);
						}}
					>
						<ListItemIcon>
							<ManageAccountsIcon fontSize="small" />
						</ListItemIcon>
						<ListItemText>Manage Account</ListItemText>
					</MenuItem>

					<MenuItem
						onClick={() => {
							handleProfileClose();
							useAppStore.getState().resetLayout();
							const api = useDockviewApi.getState().api;
							if (api) resetDockview(api);
						}}
					>
						<ListItemText>Reset layout</ListItemText>
					</MenuItem>

					{/* Admin-specific menu item */}
					{isAdmin && (
						<MenuItem
							onClick={() => {
								handleProfileClose();
								setAdminPanelOpen(true);
							}}
						>
							<ListItemIcon>
								<AdminPanelSettingsIcon fontSize="small" />
							</ListItemIcon>
							<ListItemText>Admin Panel</ListItemText>
						</MenuItem>
					)}

					{/* Logout - available to all */}
					<MenuItem onClick={handleLogout}>
						<ListItemIcon>
							<LogoutIcon fontSize="small" />
						</ListItemIcon>
						<ListItemText>Logout</ListItemText>
					</MenuItem>
				</Menu>

				{/* Hidden file input for button upload */}
				<input
					type="file"
					ref={fileInputRef}
					style={{ display: "none" }}
					onChange={handleFileInputChange}
					accept={TRAJECTORY_ACCEPT}
				/>

				{/* Main content row with activity bars, sidebars, and dockview */}
				<Box sx={{ display: "flex", flexGrow: 1, minHeight: 0 }}>
					<ActivityBar position="left" />
					<SidebarZone position="left" />
					<DockviewLayout />
					<SidebarZone position="right" />
					<ActivityBar position="right" />
				</Box>

				<BottomZone />
				<ActivityBar position="bottom" />

				<FrameProgressBar />
			</Box>

			<ConnectionDialog
				open={connectionDialogOpen}
				onClose={() => setConnectionDialogOpen(false)}
			/>

			<LoginDialog
				open={loginDialogOpen}
				onClose={() => setLoginDialogOpen(false)}
			/>

			<RegisterDialog
				open={registerDialogOpen}
				onClose={() => setRegisterDialogOpen(false)}
			/>

			<AdminPanel
				open={adminPanelOpen}
				onClose={() => setAdminPanelOpen(false)}
			/>

			<UserProfileDialog
				open={userProfileDialogOpen}
				onClose={() => setUserProfileDialogOpen(false)}
			/>

			<SiMGenTutorialDialog
				open={tutorialDialogOpen}
				onClose={() => setTutorialDialogOpen(false)}
				url="https://slides.com/rokasel/zndrawtutorial-9cc179/fullscreen?style=light"
			/>

			{snackbar && (
				<Snackbar
					open={snackbar.open}
					autoHideDuration={6000}
					onClose={hideSnackbar}
					anchorOrigin={{ vertical: "top", horizontal: "center" }}
				>
					<Alert severity={snackbar.severity} onClose={hideSnackbar}>
						{snackbar.message}
					</Alert>
				</Snackbar>
			)}

			{/* Progress notifications */}
			<ProgressNotifications />
		</>
	);
}
