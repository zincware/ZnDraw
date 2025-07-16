import { Box, Button } from "@mui/material";
import { useEffect, useRef, useState } from "react";
import { FaFileCirclePlus, FaFilm, FaRocket } from "react-icons/fa6";
import { TbPlugConnected } from "react-icons/tb";
import * as znsocket from "znsocket";
import { client, socket } from "../../socket";
import { BtnTooltip } from "../tooltips";

export function SiMGenButtons({
	visible,
	token,
	tutorialURL,
	setTutorialModalShow,
}: {
	visible: boolean;
	token: string;
	tutorialURL?: string;
	setTutorialModalShow: (show: boolean) => void;
}) {
	const [disabledBtn, setDisabledBtn] = useState<boolean>(false);
	const queueRef = useRef<any>(null);

	useEffect(() => {
		const queue = new znsocket.Dict({
			client: client,
			key: `queue:${token}:modifier`,
		});
		queueRef.current = queue;

		queue.length().then((length: any) => {
			setDisabledBtn(length > 0);
		});
		queue.onRefresh(async (x: any) => {
			const length = await queue.length();
			setDisabledBtn(length > 0);
		});

		return () => {
			queue.offRefresh();
		};
	}, [token]);

	const runConnect = () => {
		if (queueRef.current) {
			queueRef.current.Connect = {};
			socket.emit("room:worker:run");
			setDisabledBtn(true);
		}
	};

	const runGenerate = () => {
		if (queueRef.current) {
			queueRef.current.SiMGenDemo = {};
			socket.emit("room:worker:run");
			setDisabledBtn(true);
		}
	};

	const createNewCanvas = () => {
		if (queueRef.current) {
			queueRef.current.NewCanvas = {};
			socket.emit("room:worker:run");
			setDisabledBtn(true);
		}
	};

	return (
		<>
			{visible && (
				<Box
					sx={{
						display: "flex",
						alignItems: "center",
						gap: 0.5,
						border: 1.5,
						borderColor: "success.main",
						borderRadius: 1,
						px: 2,
						py: 0.5,
						mr: 1,
						backgroundColor: "transparent",
					}}
				>
					<Button
						href="https://github.com/RokasEl/simgen"
						target="_blank"
						sx={{ 
							textTransform: "none", 
							color: "success.main", 
							minWidth: "auto",
							p: 0.25,
							fontSize: "1rem",
							fontWeight: 600
						}}
					>
						SiMGen
					</Button>
					<BtnTooltip text="Connect selected atoms (shift click)">
						<Button
							variant="contained"
							color="success"
							onClick={runConnect}
							disabled={disabledBtn}
							size="small"
							sx={{ 
								mr: 0.5,
								minWidth: 'auto',
								px: 1,
								gap: 0.5,
								textTransform: 'none'
							}}
						>
							<TbPlugConnected /> Connect
						</Button>
					</BtnTooltip>
					<BtnTooltip text="Run SiMGen molecular generation">
						<Button
							variant="contained"
							color="success"
							onClick={runGenerate}
							disabled={disabledBtn}
							size="small"
							sx={{ 
								mr: 0.5,
								minWidth: 'auto',
								px: 1,
								gap: 0.5,
								textTransform: 'none'
							}}
						>
							<FaRocket /> Generate
						</Button>
					</BtnTooltip>
					<BtnTooltip text="Replace scene with empty canvas and enter drawing mode">
						<Button
							variant="contained"
							color="success"
							onClick={createNewCanvas}
							disabled={disabledBtn}
							size="small"
							sx={{ 
								mr: 0.5,
								minWidth: 'auto',
								px: 1,
								gap: 0.5,
								textTransform: 'none'
							}}
						>
							<FaFileCirclePlus /> New Canvas
						</Button>
					</BtnTooltip>
					{tutorialURL && (
						<BtnTooltip text="Tutorial">
							<Button
								variant="outlined"
								color="success"
								onClick={() => setTutorialModalShow(true)}
								startIcon={<FaFilm />}
								size="small"
								sx={{ 
									minWidth: 'auto',
									px: 1,
									textTransform: 'none',
								}}
							>
								Tutorial
							</Button>
						</BtnTooltip>
					)}
				</Box>
			)}
		</>
	);
}
