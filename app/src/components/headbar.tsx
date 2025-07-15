import { AppBar, Box, Button, IconButton, ToggleButton, Toolbar, Typography } from "@mui/material";
import { useEffect, useMemo, useState } from "react";
import { version } from "../../package.json";
import "katex/dist/katex.min.css";
import {
	FaCloudDownloadAlt,
	FaCode,
	FaDownload,
	FaFilm,
	FaHandSparkles,
	FaLock,
	FaLockOpen,
	FaMoon,
	FaSun,
	FaTerminal,
} from "react-icons/fa";
import { GrHelpBook } from "react-icons/gr";
import { socket } from "../socket";
import { BtnTooltip } from "./tooltips";

import { FaArrowRotateRight, FaPencil } from "react-icons/fa6";
import { MdAddChart, MdExitToApp } from "react-icons/md";
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
	
	const iconSize = "1.4rem";

	return (
		<>
			<AppBar position="fixed" sx={{ height: 50, bgcolor: "background.paper" }}>
				<Toolbar sx={{ height: 50, minHeight: 50, justifyContent: "space-between" }}>
					{/* --- LEFT ALIGNED ITEMS --- */}
					<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
						<Button
							href="https://github.com/zincware/zndraw"
							target="_blank"
							sx={{ textTransform: "none", color: "text.primary" }}
						>
							<Typography variant="h6" component="div">
								ZnDraw {version}
							</Typography>
						</Button>
						{showSiMGen && (
							<>
								<Typography variant="h6" sx={{ color: "text.secondary" }}>+</Typography>
								<Button
									href="https://github.com/RokasEl/simgen"
									target="_blank"
									sx={{ textTransform: "none", color: "text.primary" }}
								>
									<Typography variant="h6" component="div">
										SiMGen
									</Typography>
								</Button>
							</>
						)}
						{/* Spacer */}
						<Box sx={{ width: 24 }} />

						<BtnTooltip text="Reset Scene">
							<IconButton color="error" onClick={() => setRefreshModalShow(true)}>
								<FaArrowRotateRight style={{ fontSize: iconSize }} />
							</IconButton>
						</BtnTooltip>
						<BtnTooltip text="Activate Drawing Tool">
							<ToggleButton
								value="1"
								selected={isDrawing}
								onChange={() => setIsDrawing((prev: boolean) => !prev)}
								size="small"
								sx={{ border: 0 }} // Unify style with IconButton
							>
								<FaPencil style={{ fontSize: iconSize }} />
							</ToggleButton>
						</BtnTooltip>
						<BtnTooltip text="Remove all guiding points and geometries">
							<IconButton color="primary" onClick={handleRemovePointsGeometries}>
								<FaHandSparkles style={{ fontSize: iconSize }} />
							</IconButton>
						</BtnTooltip>
					</Box>

					{/* --- RIGHT ALIGNED ITEMS --- */}
					<Box sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
						<SiMGenButtons visible={showSiMGen} token={token} />

						{tutorialURL && (
							<BtnTooltip text="Tutorial">
								<Button
									variant="text" // More consistent with icon buttons
									color="warning"
									onClick={() => setTutorialModalShow(true)}
									startIcon={<FaFilm />}
								>
									Tutorial
								</Button>
							</BtnTooltip>
						)}
						<BtnTooltip text="Access console">
							<ToggleButton
								value="1"
								selected={consoleShow}
								onChange={() => setConsoleShow((prev: boolean) => !prev)}
								size="small"
								sx={{ border: 0 }} // Unify style with IconButton
							>
								<FaTerminal style={{ fontSize: iconSize }} />
							</ToggleButton>
						</BtnTooltip>
						<BtnTooltip text="Python access">
							<IconButton color="primary" onClick={() => setConnectModalShow(true)}>
								<FaCode style={{ fontSize: iconSize }} />
							</IconButton>
						</BtnTooltip>
						<BtnTooltip text="Add plots window">
							<IconButton
								color="primary"
								onClick={() => setAddPlotsWindow((prev: number) => prev + 1)}
							>
								<MdAddChart style={{ fontSize: iconSize }} />
							</IconButton>
						</BtnTooltip>
						<FileUpload />
						<BtnTooltip text="Download">
							<IconButton
								color="primary"
								href={`${basePath}download`}
								target="_blank"
							>
								<FaDownload style={{ fontSize: iconSize }} />
							</IconButton>
						</BtnTooltip>
						{zntrackAvailable && (
							<BtnTooltip text="Open File via DVC">
								<IconButton color="primary" onClick={() => setRemoteFileModalShow(true)}>
									<FaCloudDownloadAlt style={{ fontSize: iconSize }} />
								</IconButton>
							</BtnTooltip>
						)}
						<BtnTooltip text="Help">
							<IconButton color="primary" onClick={() => setHelpModalShow(true)}>
								<GrHelpBook style={{ fontSize: iconSize }} />
							</IconButton>
						</BtnTooltip>
						<BtnTooltip text="Switch Colormode">
							<IconButton color="primary" onClick={handleColorMode}>
								{colorMode === "light" ? (
									<FaSun style={{ fontSize: iconSize }} />
								) : (
									<FaMoon style={{ fontSize: "1.2rem" }} />
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
					</Box>
				</Toolbar>
			</AppBar>
            
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
