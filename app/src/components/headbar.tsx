import {
	AppBar,
	Box,
	Button,
	IconButton,
	ToggleButton,
	Toolbar,
	Typography,
} from "@mui/material";
import "katex/dist/katex.min.css";
import { useEffect, useMemo, useState } from "react";
import {
	FaCloudDownloadAlt,
	FaCode,
	FaDownload,
	FaHandSparkles,
	FaLock,
	FaLockOpen,
	FaMoon,
	FaSun,
	FaTerminal,
	FaUpload,
} from "react-icons/fa";
import { FaArrowRotateRight, FaPencil } from "react-icons/fa6";
import { GrHelpBook } from "react-icons/gr";
import { MdAddChart, MdExitToApp } from "react-icons/md";
import { version } from "../../package.json";
import { socket } from "../socket";
import { BtnTooltip } from "./tooltips";
import { ConnectModal } from "./headbar/ConnectModal";
import { ConsoleWindow } from "./headbar/ConsoleWindow";
import { FileUpload } from "./headbar/FileUpload";
import { HelpModel } from "./headbar/HelpModal";
import { RefreshModal } from "./headbar/RefreshModal";
import { RemoteFileModal } from "./headbar/RemoteFileModal";
import { SiMGenButtons } from "./headbar/SiMGenButtons";
import { TutorialModal } from "./headbar/TutorialModal";

interface HeadBarProps {
	room: string;
	colorMode: string;
	handleColorMode: any;
	setIsDrawing: any;
	setGeometries: any;
	setPoints: any;
	isDrawing: boolean;
	tutorialURL: string;
	showSiMGen: boolean;
	modifierQueue: number;
	isAuthenticated: boolean;
	roomLock: boolean;
	setAddPlotsWindow: any;
	messages: any[];
	setMessages: any;
	token: string;
	step: number;
	selection: Set<number>;
}

const HeadBar = ({
	room,
	colorMode,
	handleColorMode,
	setIsDrawing,
	setGeometries,
	setPoints,
	isDrawing,
	tutorialURL,
	showSiMGen,
	modifierQueue,
	isAuthenticated,
	roomLock,
	setAddPlotsWindow,
	messages,
	setMessages,
	token,
	step,
	selection,
}: HeadBarProps) => {
	const [helpModalShow, setHelpModalShow] = useState(false);
	const [connectModalShow, setConnectModalShow] = useState(false);
	const [refreshModalShow, setRefreshModalShow] = useState(false);
	const [tutorialModalShow, setTutorialModalShow] = useState(false);
	const [consoleShow, setConsoleShow] = useState(false);
	const [remoteFileModalShow, setRemoteFileModalShow] = useState(false);
	const [zntrackAvailable, setZntrackAvailable] = useState(false);

	useEffect(() => {
		socket.emit("zntrack:available", (available: boolean) => {
			setZntrackAvailable(available);
		});
	}, []);

	useEffect(() => {
		setConsoleShow(showSiMGen);
	}, [showSiMGen]);

	const handleRemovePointsGeometries = () => {
		socket.emit("room:geometry:set", []);
		socket.emit("room:point:set", []);
		setGeometries([]);
		setPoints([]);
	};
	const basePath = useMemo(() => import.meta.env.BASE_URL || "/", []);
	const iconSize = "1.25rem";

	// Common styles for icon buttons to ensure consistent alignment
	const iconButtonStyle = {
		width: 40,
		height: 40,
		display: "flex",
		alignItems: "center",
		justifyContent: "center",
	};

	const HeadBarHeight = 50; // Height of the AppBar

	return (
		<>
			<AppBar
				position="fixed"
				sx={{
					height: HeadBarHeight,
					bgcolor: "background.paper",
					top: 0,
					left: 0,
					right: 0,
					zIndex: 1100,
				}}
			>
				<Toolbar
					sx={{
						height: HeadBarHeight,
						minHeight: HeadBarHeight,
						display: "flex",
						alignItems: "center",
						gap: 0.5,
						// Override Material-UI's default toolbar behavior
						"&.MuiToolbar-root": {
							minHeight: HeadBarHeight,
							height: HeadBarHeight,
							padding: 0,
							paddingLeft: 2,
							paddingRight: 2,
						},
					}}
				>
					{/* --- LEFT ALIGNED ITEMS --- */}
					<Button
						href="https://github.com/zincware/zndraw"
						target="_blank"
						sx={{ textTransform: "none", color: "text.primary", mr: 1 }}
					>
						<Typography variant="h6" component="div">
							ZnDraw {version}
						</Typography>
					</Button>

					<BtnTooltip text="Reset Scene">
						<IconButton
							color="error"
							onClick={() => setRefreshModalShow(true)}
							sx={iconButtonStyle}
						>
							<FaArrowRotateRight style={{ fontSize: iconSize }} />
						</IconButton>
					</BtnTooltip>
					<BtnTooltip text="Activate Drawing Tool">
						<ToggleButton
							value="1"
							selected={isDrawing}
							onChange={() => setIsDrawing((prev: boolean) => !prev)}
							size="small"
							sx={{
								border: 0,
								...iconButtonStyle,
								"&.MuiToggleButton-root": {
									"&:hover": {
										backgroundColor: (theme) => theme.palette.action.hover,
									},
								},
							}}
						>
							<FaPencil style={{ fontSize: iconSize }} />
						</ToggleButton>
					</BtnTooltip>
					<BtnTooltip text="Remove all guiding points and geometries">
						<IconButton
							color="primary"
							onClick={handleRemovePointsGeometries}
							sx={iconButtonStyle}
						>
							<FaHandSparkles style={{ fontSize: iconSize }} />
						</IconButton>
					</BtnTooltip>

					{/* This spacer pushes all subsequent items to the right */}
					<Box sx={{ flexGrow: 1 }} />

					{/* --- RIGHT ALIGNED ITEMS --- */}
					<SiMGenButtons 
						visible={showSiMGen} 
						token={token} 
						tutorialURL={tutorialURL}
						setTutorialModalShow={setTutorialModalShow}
					/>
					<BtnTooltip text="Access console">
						<ToggleButton
							value="1"
							selected={consoleShow}
							onChange={() => setConsoleShow((prev: boolean) => !prev)}
							size="small"
							sx={{
								border: 0,
								...iconButtonStyle,
								"&.MuiToggleButton-root": {
									"&:hover": {
										backgroundColor: (theme) => theme.palette.action.hover,
									},
								},
							}}
						>
							<FaTerminal style={{ fontSize: iconSize }} />
						</ToggleButton>
					</BtnTooltip>
					<BtnTooltip text="Python access">
						<IconButton
							color="primary"
							onClick={() => setConnectModalShow(true)}
							sx={iconButtonStyle}
						>
							<FaCode style={{ fontSize: iconSize }} />
						</IconButton>
					</BtnTooltip>
					<BtnTooltip text="Add plots window">
						<IconButton
							color="primary"
							onClick={() => setAddPlotsWindow((prev: number) => prev + 1)}
							sx={iconButtonStyle}
						>
							<MdAddChart style={{ fontSize: iconSize }} />
						</IconButton>
					</BtnTooltip>
					<FileUpload
						renderButton={(handleClick) => (
							<BtnTooltip text="Upload File">
								<IconButton
									color="primary"
									onClick={handleClick}
									component="span"
									sx={iconButtonStyle}
								>
									<FaUpload style={{ fontSize: iconSize }} />
								</IconButton>
							</BtnTooltip>
						)}
					/>
					<BtnTooltip text="Download">
						<IconButton
							color="primary"
							href={`${basePath}download`}
							target="_blank"
							sx={iconButtonStyle}
						>
							<FaDownload style={{ fontSize: iconSize }} />
						</IconButton>
					</BtnTooltip>
					{zntrackAvailable && (
						<BtnTooltip text="Open File via DVC">
							<IconButton
								color="primary"
								onClick={() => setRemoteFileModalShow(true)}
							>
								<FaCloudDownloadAlt style={{ fontSize: iconSize }} />
							</IconButton>
						</BtnTooltip>
					)}
					<BtnTooltip text="Help">
						<IconButton
							color="primary"
							onClick={() => setHelpModalShow(true)}
							sx={iconButtonStyle}
						>
							<GrHelpBook style={{ fontSize: iconSize }} />
						</IconButton>
					</BtnTooltip>
					<BtnTooltip text="Switch Colormode">
						<IconButton
							color="primary"
							onClick={handleColorMode}
							sx={iconButtonStyle}
						>
							{colorMode === "light" ? (
								<FaSun style={{ fontSize: iconSize }} />
							) : (
								<FaMoon style={{ fontSize: iconSize }} />
							)}
						</IconButton>
					</BtnTooltip>
					{isAuthenticated && (
						<>
							<BtnTooltip
								text={roomLock ? "Unlock this room" : "Lock this room"}
							>
								<IconButton
									color="error"
									onClick={() => socket.emit("room:lock:set", !roomLock)}
								>
									{roomLock ? (
										<FaLock style={{ fontSize: iconSize }} />
									) : (
										<FaLockOpen style={{ fontSize: iconSize }} />
									)}
								</IconButton>
							</BtnTooltip>
							<BtnTooltip text="Close ZnDraw">
								<IconButton
									color="error"
									onClick={() => {
										socket.emit("shutdown");
										close();
									}}
								>
									<MdExitToApp style={{ fontSize: iconSize }} />
								</IconButton>
							</BtnTooltip>
						</>
					)}
				</Toolbar>
			</AppBar>
			{/* --- Modals and other components --- */}
			<HelpModel show={helpModalShow} onHide={() => setHelpModalShow(false)} />
			<ConnectModal
				show={connectModalShow}
				onHide={() => setConnectModalShow(false)}
				token={token}
			/>
			<RefreshModal
				show={refreshModalShow}
				onHide={() => setRefreshModalShow(false)}
				room={room}
				token={token}
			/>
			<TutorialModal
				show={tutorialModalShow}
				onHide={() => setTutorialModalShow(false)}
				url={tutorialURL}
			/>
			<RemoteFileModal
				show={remoteFileModalShow}
				onHide={() => setRemoteFileModalShow(false)}
				colorMode={colorMode}
			/>
			{consoleShow && (
				<ConsoleWindow
					messages={messages}
					setConsoleShow={setConsoleShow}
					token={token}
					setMessages={setMessages}
					colorMode={colorMode}
					step={step}
					selection={selection}
				/>
			)}
		</>
	);
};

export default HeadBar;
