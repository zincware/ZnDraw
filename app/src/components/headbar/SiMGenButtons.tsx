import { Button } from "@mui/material";
import { useEffect, useRef, useState } from "react";
import { FaFileCirclePlus, FaRocket } from "react-icons/fa6";
import { TbPlugConnected } from "react-icons/tb";
import * as znsocket from "znsocket";
import { client, socket } from "../../socket";
import { BtnTooltip } from "../tooltips";

export function SiMGenButtons({
	visible,
	token,
}: {
	visible: boolean;
	token: string;
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
				<>
					<BtnTooltip text="Connect selected atoms (shift click)">
						<Button
							variant="success"
							className="mx-1"
							onClick={runConnect}
							disabled={disabledBtn}
						>
							<TbPlugConnected /> Connect
						</Button>
					</BtnTooltip>
					<BtnTooltip text="Run SiMGen molecular generation">
						<Button
							variant="success"
							className="mx-1"
							onClick={runGenerate}
							disabled={disabledBtn}
						>
							<FaRocket /> Generate
						</Button>
					</BtnTooltip>
					<BtnTooltip text="Replace scene with empty canvas and enter drawing mode">
						<Button
							variant="success"
							className="mx-1"
							onClick={createNewCanvas}
							disabled={disabledBtn}
						>
							<FaFileCirclePlus /> New Canvas
						</Button>
					</BtnTooltip>
				</>
			)}
		</>
	);
}
