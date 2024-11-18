import { useRef, useState } from "react";
import { useEffect } from "react";
import { Button, ButtonGroup, Card, Form, Nav, Navbar } from "react-bootstrap";
import {
	FaCube,
	FaGithub,
	FaPlay,
	FaRegChartBar,
	FaRegHandPointer,
	FaRegMap,
} from "react-icons/fa";
import { FaCircleNodes, FaGear } from "react-icons/fa6";
import { IoStop } from "react-icons/io5";
import Select from "react-select";
import * as znsocket from "znsocket";
import { client, socket } from "../socket";
import { BtnTooltip } from "./tooltips";

import { JSONEditor } from "@json-editor/json-editor";

JSONEditor.defaults.options.theme = "bootstrap5";
JSONEditor.defaults.options.iconlib = "fontawesome5";
JSONEditor.defaults.options.object_background = "bg-body-color";
JSONEditor.defaults.options.disable_edit_json = true;
JSONEditor.defaults.options.disable_properties = true;
JSONEditor.defaults.options.disable_collapse = true;
JSONEditor.defaults.options.no_additional_properties = true;
JSONEditor.defaults.options.keep_oneof_values = false;
JSONEditor.defaults.editors.object.options.titleHidden = true;

interface SidebarMenuProps {
	visible: boolean;
	closeMenu: () => void;
	token: string;
	name: string;
	sendImmediately: boolean;
}

const SidebarMenu = ({
	visible,
	closeMenu,
	token,
	name,
	sendImmediately,
}: SidebarMenuProps) => {
	const [userInput, setUserInput] = useState<string>(undefined);
	const [schema, setSchema] = useState<any>({});
	const [sharedSchema, setSharedSchema] = useState<any>({});
	const [editorValue, setEditorValue] = useState<any>(null);
	const [disabledBtn, setDisabledBtn] = useState<boolean>(false);
	const initialTrigger = useRef<boolean>(true);
	const editorRef = useRef<HTMLDivElement>(null);
	const queueRef = useRef<any>(null);

	useEffect(() => {
		if (visible) {
			socket.emit("schema:refresh");
			initialTrigger.current = true;
		}
	}, [visible]);

	useEffect(() => {
		const con = new znsocket.Dict({
			client: client,
			key: `schema:${token}:${name}`,
		});

		const sharedCon = new znsocket.Dict({
			client: client,
			key: `schema:default:${name}`,
		});

		const queue = new znsocket.Dict({
			client: client,
			key: `queue:${token}:${name}`,
		});
		queueRef.current = queue;

		// initial load
		con.entries().then((items: any) => {
			const result = Object.fromEntries(items);
			setSchema(result);
		});
		sharedCon.entries().then((items: any) => {
			const result = Object.fromEntries(items);
			setSharedSchema(result);
		});

		queue.length().then((length: any) => {
			setDisabledBtn(length > 0);
		});
		queue.onRefresh(async (x: any) => {
			const length = await queue.length();
			setDisabledBtn(length > 0);
		});

		con.onRefresh(async (x: any) => {
			const items = await con.entries();
			const result = Object.fromEntries(items);
			setSchema(result);
		});

		sharedCon.onRefresh(async (x: any) => {
			const items = await sharedCon.entries();
			const result = Object.fromEntries(items);
			setSharedSchema(result);
		});

		return () => {
			con.offRefresh();
			sharedCon.offRefresh();
			queue.offRefresh();
		};
	}, [token]);

	// set the default userInput to the first key in the schema, if userInput is empty
	useEffect(() => {
		if (
			userInput === undefined &&
			Object.keys({ ...sharedSchema, ...schema }).length > 0
		) {
			setUserInput(Object.keys({ ...sharedSchema, ...schema })[0]);
		}
	}, [schema, sharedSchema, userInput]);

	useEffect(() => {
		let editor: any;
		const fullSchema = { ...sharedSchema, ...schema };
		if (editorRef.current && fullSchema[userInput]) {
			editor = new JSONEditor(editorRef.current, {
				schema: fullSchema[userInput],
			});

			editor.on("ready", () => {
				if (editor.validate()) {
					const editorValue = editor.getValue();
					setEditorValue(editorValue);
				}
			});

			editor.on("change", () => {
				if (editor.ready) {
					if (editor.validate()) {
						const editorValue = editor.getValue();
						setEditorValue(editorValue);
					}
				}
			});
		}
		return () => {
			if (editor) {
				editor.destroy();
				setEditorValue(null);
			}
		};
	}, [userInput, schema, sharedSchema]);

	function submitEditor() {
		if (editorValue && userInput && queueRef.current) {
			setDisabledBtn(true);
			// queueRef.current.push({ [userInput]: editorValue });
			queueRef.current[userInput] = editorValue;
			socket.emit("room:worker:run");
		}
	}

	useEffect(() => {
		if (sendImmediately && editorValue && userInput && queueRef.current) {
			if (initialTrigger.current) {
				initialTrigger.current = false;
				return;
			}
			submitEditor();
		}
	}, [editorValue]);

	function cancelTask() {
		if (queueRef.current) {
			queueRef.current.clear().then(() => {
				setDisabledBtn(false);
			});
		}
	}

	function capitalizeFirstLetter(text: string) {
		if (!text) return ""; // Handle empty or undefined text
		return text.charAt(0).toUpperCase() + text.slice(1);
	}

	return (
		<Card
			className="rounded-0 border-start-0 overflow-y-auto rounded-end"
			style={{
				position: "fixed",
				top: "50px",
				left: "50px",
				height: "100%",
				maxWidth: "35%",
				margin: 0,
				padding: 0,
				display: visible ? "block" : "none",
			}}
		>
			<Card.Header
				className="d-flex justify-content-between align-items-center"
				style={{
					backgroundColor: "inherit", // Use the same background color as the rest of the card
				}}
			>
				<Card.Title>{capitalizeFirstLetter(name)}</Card.Title>
				<Button variant="close" className="ms-auto" onClick={closeMenu} />
			</Card.Header>
			<Card.Body style={{ paddingBottom: 80 }}>
				<Form.Group className="d-flex align-items-center">
					<Form.Select
						aria-label="Default select example"
						onChange={(e) => setUserInput(e.target.value)}
						value={userInput}
					>
						<option />
						{Object.keys({ ...sharedSchema, ...schema }).map((key) => (
							<option key={key} value={key}>
								{key}
							</option>
						))}
					</Form.Select>
					{!sendImmediately && (
						<ButtonGroup aria-label="Basic example">
							<Button
								variant="outline-primary"
								onClick={submitEditor}
								className="ms-2 d-flex align-items-center"
								disabled={disabledBtn}
							>
								<FaPlay className="me-1" /> Submit
							</Button>
							<Button variant="outline-danger" onClick={cancelTask}>
								<IoStop />
							</Button>
						</ButtonGroup>
					)}
				</Form.Group>
				<div ref={editorRef} />
			</Card.Body>
		</Card>
	);
};

function SideBar({ token }: { token: string }) {
	const [visibleOption, setVisibleOption] = useState<string>("");
	useEffect(() => {
		if (visibleOption !== "") {
			// emit visibleOption:schema e.g. selection:schema
			socket.emit(`${visibleOption}:schema`);
		}
	}, [visibleOption]);

	// if any menu is open and you click escape, close it
	useEffect(() => {
		function handleKeyDown(event: KeyboardEvent) {
			if (event.key === "Escape") {
				setVisibleOption("");
			}
		}

		document.addEventListener("keydown", handleKeyDown);

		return () => {
			document.removeEventListener("keydown", handleKeyDown);
		};
	}, []);

	return (
		<>
			<Navbar className="flex-column bg-body-primary bg-tertiary">
				<Card
					border="tertiary py-0 px-0"
					className="rounded-0 border-top-0"
					style={{
						position: "fixed",
						top: "50px",
						left: "0",
						height: "100%",
						width: "50px",
					}}
				>
					<BtnTooltip text="Selection" placement="right">
						<Nav className="mx-auto my-1">
							<Button
								variant="outline-tertiary"
								onClick={() =>
									setVisibleOption(
										visibleOption !== "selection" ? "selection" : "",
									)
								}
							>
								<FaRegHandPointer />
							</Button>
						</Nav>
					</BtnTooltip>
					<BtnTooltip text="Interaction" placement="right">
						<Nav className="mx-auto my-1">
							<Button
								variant="outline-tertiary"
								onClick={() =>
									setVisibleOption(
										visibleOption !== "modifier" ? "modifier" : "",
									)
								}
							>
								<FaCircleNodes />
							</Button>
						</Nav>
					</BtnTooltip>
					<BtnTooltip text="Settings" placement="right">
						<Nav className="mx-auto my-1">
							<Button
								variant="outline-tertiary"
								onClick={() =>
									setVisibleOption(
										visibleOption !== "settings" ? "settings" : "",
									)
								}
							>
								<FaGear />
							</Button>
						</Nav>
					</BtnTooltip>
					<BtnTooltip text="Geometry" placement="right">
						<Nav className="mx-auto my-1">
							<Button
								variant="outline-tertiary"
								onClick={() =>
									setVisibleOption(
										visibleOption !== "geometry" ? "geometry" : "",
									)
								}
							>
								<FaRegMap />
							</Button>
						</Nav>
					</BtnTooltip>
					<BtnTooltip text="Analysis" placement="right">
						<Nav className="mx-auto my-1">
							<Button
								variant="outline-tertiary"
								onClick={() =>
									setVisibleOption(
										visibleOption !== "analysis" ? "analysis" : "",
									)
								}
							>
								<FaRegChartBar />
							</Button>
						</Nav>
					</BtnTooltip>
					<BtnTooltip text="View on GitHub" placement="right">
						<Nav className="mx-auto my-1">
							<Button
								variant="outline-tertiary"
								href="https://github.com/zincware/ZnDraw"
								target="_blank"
							>
								<FaGithub />
							</Button>
						</Nav>
					</BtnTooltip>
				</Card>
			</Navbar>
			<SidebarMenu
				name="selection"
				visible={visibleOption === "selection"} // remove
				token={token}
				closeMenu={() => setVisibleOption("")}
				sendImmediately={false}
			/>
			<SidebarMenu
				name="modifier"
				visible={visibleOption === "modifier"}
				token={token}
				closeMenu={() => setVisibleOption("")}
				sendImmediately={false}
			/>
			<SidebarMenu
				name="settings"
				visible={visibleOption === "settings"}
				token={token}
				closeMenu={() => setVisibleOption("")}
				sendImmediately={true}
			/>
			<SidebarMenu
				name="geometry"
				visible={visibleOption === "geometry"}
				token={token}
				closeMenu={() => setVisibleOption("")}
				sendImmediately={false}
			/>
			<SidebarMenu
				name="analysis"
				visible={visibleOption === "analysis"}
				token={token}
				closeMenu={() => setVisibleOption("")}
				sendImmediately={false}
			/>
		</>
	);
}

export default SideBar;
