import { AppBar, Box, Button, ToggleButton, Toolbar } from "@mui/material";
import { useEffect, useMemo, useState } from "react";
import { version } from "../../package.json";
import "katex/dist/katex.min.css"; // `rehype-katex` does not import the CSS for you
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
		console.log("remove points and geometries");
		socket.emit("room:geometry:set", []);
		socket.emit("room:point:set", []);
		setGeometries([]);
		setPoints([]);
	};
	const basePath = useMemo(() => import.meta.env.BASE_URL || "/", []);

	return (
		<>
			<AppBar position="fixed" sx={{ height: 50, bgcolor: "background.paper" }}>
				<Toolbar>
					<Box sx={{ flexGrow: 1 }}>
						<Button
							size="large"
							color="inherit"
							href="https://github.com/zincware/zndraw"
							target="_blank"
						>
							ZnDraw {version}
						</Button>
						{showSiMGen && (
							<>
								<Button size="large" color="inherit">
									+
								</Button>
								<Button
									size="large"
									color="inherit"
									href="https://github.com/RokasEl/simgen"
									target="_blank"
								>
									SiMGen
								</Button>
							</>
						)}
					</Box>
					<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
						<BtnTooltip text="Reset Scene">
							<Button
								variant="outlined"
								color="error"
								size="small"
								onClick={() => {
									setRefreshModalShow((prev: boolean) => !prev);
								}}
							>
								<FaArrowRotateRight />
							</Button>
						</BtnTooltip>
						<BtnTooltip text="Activate Drawing Tool">
							<ToggleButton
								value="1"
								selected={isDrawing}
								onChange={(e) => {
									setIsDrawing((prev: boolean) => !prev);
								}}
							>
								<FaPencil />
							</ToggleButton>
						</BtnTooltip>
						<BtnTooltip text="Remove all guiding points and geometries">
							<Button
								variant="outlined"
								color="primary"
								size="small"
								onClick={handleRemovePointsGeometries}
							>
								<FaHandSparkles />
							</Button>
						</BtnTooltip>
						<SiMGenButtons visible={showSiMGen} token={token} />
					</Box>
					<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
						{tutorialURL && (
							<BtnTooltip text="Tutorial">
								<Button
									variant="warning"
									className="mx-1"
									onClick={() => setTutorialModalShow(true)}
								>
									Tutorial <FaFilm />
								</Button>
							</BtnTooltip>
						)}
						<BtnTooltip text="Access console">
							<ToggleButton
								value="1"
								selected={consoleShow}
								onChange={() => {
									setConsoleShow((prev: boolean) => !prev);
								}}
							>
								<FaTerminal />
							</ToggleButton>
						</BtnTooltip>
						<BtnTooltip text="Python access">
							<Button
								variant="outline-primary"
								className="mx-1"
								onClick={() => setConnectModalShow(true)}
							>
								<FaCode />
							</Button>
						</BtnTooltip>
						<BtnTooltip text="Add plots window">
							<Button
								variant="outline-primary"
								className="mx-1"
								onClick={() => setAddPlotsWindow((prev: number) => prev + 1)}
							>
								<MdAddChart />
							</Button>
						</BtnTooltip>
						<FileUpload />
						<BtnTooltip text="Download">
							<Button
								variant="outline-primary"
								className="mx-1"
								href={`${basePath}download`}
								target="_blank"
							>
								<FaDownload />
							</Button>
						</BtnTooltip>
						{zntrackAvailable && (
							<BtnTooltip text="Open File via DVC">
								<Button
									variant="outline-primary"
									className="mx-1"
									onClick={() => {
										setRemoteFileModalShow(true);
									}}
								>
									<FaCloudDownloadAlt />
								</Button>
							</BtnTooltip>
						)}
						<BtnTooltip text="Help">
							<Button
								variant="outline-primary"
								className="mx-1"
								onClick={() => setHelpModalShow(true)}
							>
								<GrHelpBook />
							</Button>
						</BtnTooltip>
						<BtnTooltip text="Switch Colormode">
							<Button
								variant="outline-primary"
								className="mx-1"
								onClick={handleColorMode}
							>
								{colorMode === "light" ? <FaSun /> : <FaMoon />}
							</Button>
						</BtnTooltip>
						{isAuthenticated && (
							<>
								<BtnTooltip
									text={roomLock ? "Unlock this room" : "Lock this room"}
								>
									<Button
										variant="outline-danger"
										className="mx-1"
										onClick={() => {
											socket.emit("room:lock:set", !roomLock);
										}}
									>
										{roomLock ? <FaLock /> : <FaLockOpen />}
									</Button>
								</BtnTooltip>
								{/* <Button variant="outline-primary" className="mx-1">
                <FaUsers />
              </Button> */}
								<BtnTooltip text="Close ZnDraw">
									<Button
										variant="outline-danger"
										className="mx-1"
										onClick={() => {
											socket.emit("shutdown");
											close();
										}}
									>
										<MdExitToApp />
									</Button>
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
