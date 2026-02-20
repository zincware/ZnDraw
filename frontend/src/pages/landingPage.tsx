import AccountCircleIcon from "@mui/icons-material/AccountCircle";
import AdminPanelSettingsIcon from "@mui/icons-material/AdminPanelSettings";
import BrushIcon from "@mui/icons-material/Brush";
import CameraAltIcon from "@mui/icons-material/CameraAlt";
import ChatIcon from "@mui/icons-material/Chat";
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
import Badge from "@mui/material/Badge";
import Box from "@mui/material/Box";
import CircularProgress from "@mui/material/CircularProgress";
import CssBaseline from "@mui/material/CssBaseline";
import IconButton from "@mui/material/IconButton";
import ListItemIcon from "@mui/material/ListItemIcon";
import ListItemText from "@mui/material/ListItemText";
import Menu from "@mui/material/Menu";
import MenuItem from "@mui/material/MenuItem";
import Snackbar from "@mui/material/Snackbar";
import Toolbar from "@mui/material/Toolbar";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import React, { useCallback, useState, useEffect, Suspense } from "react";
import AdminPanel from "../components/AdminPanel";
import LoginDialog from "../components/LoginDialog";
import FrameProgressBar from "../components/ProgressBar";
import RegisterDialog from "../components/RegisterDialog";
import RoomManagementMenu from "../components/RoomManagementMenu";
import SiMGenButtons from "../components/SiMGenButtons";
import SiMGenTutorialDialog from "../components/SiMGenTutorialDialog";
import SideBar from "../components/SideBar";
import UserProfileDialog from "../components/UserProfileDialog";

import MyScene from "../components/Canvas";
import { useDragAndDrop } from "../hooks/useDragAndDrop";
import { useKeyboardShortcuts } from "../hooks/useKeyboardShortcuts";
import { useSocketManager } from "../hooks/useSocketManager";
// Lazy load ChatWindow - markdown/syntax highlighting only loads when chat is opened
const ChatWindow = React.lazy(() => import("../components/ChatWindow"));
import Link from "@mui/material/Link";
import { useColorScheme } from "@mui/material/styles";
import { useQueryClient } from "@tanstack/react-query";
import { useNavigate, useParams } from "react-router-dom";
import { uploadTrajectory } from "../myapi/client";
import AddPlotButton from "../components/AddPlotButton";
import ConnectionDialog from "../components/ConnectionDialog";
import DropOverlay from "../components/DropOverlay";
import ProgressNotifications from "../components/ProgressNotifications";
import WindowManager from "../components/WindowManager";
import { TRAJECTORY_ACCEPT } from "../constants/fileTypes";
import { LAYOUT_CONSTANTS } from "../constants/layout";
import { connectWithAuth } from "../socket";
import { useAppStore } from "../store";
import { logout as authLogout } from "../utils/auth";
import { downloadScreenshot } from "../utils/screenshot";

export default function MainPage() {
	const { roomId } = useParams<{ roomId: string }>();
	const setRoomId = useAppStore((state) => state.setRoomId);
	const navigate = useNavigate();

	// Set roomId in store for child components that read from it
	useEffect(() => {
		if (roomId) setRoomId(roomId);
	}, [roomId, setRoomId]);

	useSocketManager({ roomId });
	useKeyboardShortcuts();

	const chatOpen = useAppStore((state) => state.chatOpen);
	const setChatOpen = useAppStore((state) => state.setChatOpen);
	const interactionMode = useAppStore((state) => state.mode);
	const enterDrawingMode = useAppStore((state) => state.enterDrawingMode);
	const exitDrawingMode = useAppStore((state) => state.exitDrawingMode);
	const enterEditingMode = useAppStore((state) => state.enterEditingMode);
	const exitEditingMode = useAppStore((state) => state.exitEditingMode);
	const chatUnreadCount = useAppStore((state) => state.chatUnreadCount);
	const serverVersion = useAppStore((state) => state.serverVersion);
	const globalSettings = useAppStore((state) => state.globalSettings);
	const userEmail = useAppStore((state) => state.user?.email ?? null);
	const isAdmin = useAppStore((state) => state.user?.is_superuser ?? false);
	const setUser = useAppStore((state) => state.setUser);
	const showSnackbar = useAppStore((state) => state.showSnackbar);
	const snackbar = useAppStore((state) => state.snackbar);
	const hideSnackbar = useAppStore((state) => state.hideSnackbar);
	const screenshotCapture = useAppStore((state) => state.screenshotCapture);
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

	const {
		isDragging,
		handleDragOver,
		handleDragEnter,
		handleDragLeave,
		handleDrop,
	} = useDragAndDrop(handleFiles);

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
						zIndex: (theme) => theme.zIndex.drawer + 1,
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
						<Tooltip title={"Toggle chat window"}>
							<IconButton
								color="inherit"
								aria-label="toggle chat"
								onClick={() => setChatOpen(!chatOpen)}
							>
								<Badge
									badgeContent={chatUnreadCount}
									color="error"
									max={99}
									invisible={chatUnreadCount === 0 || chatOpen}
								>
									<ChatIcon />
								</Badge>
							</IconButton>
						</Tooltip>
						<AddPlotButton />
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

				{/* Main content row with sidebar and center area */}
				<Box sx={{ display: "flex", flexGrow: 1, minHeight: 0 }}>
					<SideBar />

					{/* Main Content Area with drag boundary */}
					<Box
						component="main"
						sx={{
							flexGrow: 1,
							display: "flex",
							flexDirection: "column",
							minWidth: 0,
						}}
					>
						{/* THIS IS THE CRUCIAL DRAG BOUNDARY CONTAINER */}
						<Box
							className="drag-boundary-container"
							onDragOver={handleDragOver}
							onDragEnter={handleDragEnter}
							onDragLeave={handleDragLeave}
							onDrop={handleDrop}
							sx={{
								flexGrow: 1,
								position: "relative",
								overflow: "hidden",
								display: "flex",
								flexDirection: "column",
							}}
						>
							<DropOverlay isDragging={isDragging} />
							<MyScene />
							<WindowManager />
						</Box>
					</Box>
				</Box>

				{/* Progress bar spans full width outside the drag boundary */}
				<FrameProgressBar />

				{/* Chat Window - lazy loaded */}
				{chatOpen && (
					<Suspense
						fallback={
							<Box
								sx={{
									position: "fixed",
									top: "50%",
									left: "50%",
									transform: "translate(-50%, -50%)",
									zIndex: 9999,
								}}
							>
								<CircularProgress />
							</Box>
						}
					>
						<ChatWindow open={chatOpen} onClose={() => setChatOpen(false)} />
					</Suspense>
				)}
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
