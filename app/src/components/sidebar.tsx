import {
	AccountTree,
	BarChart,
	Close,
	GitHub,
	Map,
	PanTool,
	PlayArrow,
	Settings,
	Stop,
} from "@mui/icons-material";
import {
	Alert,
	AppBar,
	Box,
	Button,
	ButtonGroup,
	Card,
	CardContent,
	CardHeader,
	FormControl,
	IconButton,
	InputLabel,
	MenuItem,
	Select,
	Toolbar,
	Tooltip,
	Typography,
} from "@mui/material";
import isEqual from "lodash.isequal";
import { useCallback, useMemo, useRef, useState } from "react";
import { useEffect } from "react";
import * as znsocket from "znsocket";
import { client, socket } from "../socket";
import JSONFormsEditor from "./JSONFormsEditor";

// Material-UI theme configuration will be handled by the app-level theme provider

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
	const [userInput, setUserInput] = useState<string | undefined>(undefined);
	const [schema, setSchema] = useState<any>({});
	const [sharedSchema, setSharedSchema] = useState<any>({});
	const [editorValue, setEditorValue] = useState<any>(null);
	const [disabledBtn, setDisabledBtn] = useState<boolean>(false);
	const [validationErrors, setValidationErrors] = useState<any[]>([]);
	const initialTrigger = useRef<boolean>(true);
	const queueRef = useRef<any>(null);
	const lastSchemaRef = useRef<any>(null);
	const debounceTimerRef = useRef<NodeJS.Timeout | null>(null);

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
		queue.onRefresh(async () => {
			const length = await queue.length();
			setDisabledBtn(length > 0);
		});

		con.onRefresh(async () => {
			const items = await con.entries();
			const result = Object.fromEntries(items);
			setSchema(result);
		});

		sharedCon.onRefresh(async () => {
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
			const keys = Object.keys({ ...sharedSchema, ...schema });
			if (keys.length > 0) {
				setUserInput(keys[0]);
			}
		}
	}, [schema, sharedSchema, userInput]);

	// Memoize the full schema to prevent unnecessary re-renders
	const fullSchema = useMemo(() => {
		const merged = { ...sharedSchema, ...schema };
		console.log("Full schema updated:", merged);
		return merged;
	}, [sharedSchema, schema]);

	// Memoize the current schema to prevent unnecessary re-renders
	const currentSchema = useMemo(() => {
		if (!userInput) return null;
		const schema = fullSchema[userInput];
		console.log("Current schema for", userInput, ":", schema);
		return schema;
	}, [fullSchema, userInput]);

	// Only recreate editor data if the schema actually changes
	useEffect(() => {
		if (currentSchema && !isEqual(currentSchema, lastSchemaRef.current)) {
			console.log("Schema changed for", userInput, "recreating editor");
			lastSchemaRef.current = currentSchema;
			// Set initial value from schema defaults if available
			if (currentSchema.properties) {
				const defaultValue: any = {};
				for (const [key, prop] of Object.entries(currentSchema.properties)) {
					const property = prop as any;
					if (property.default !== undefined) {
						defaultValue[key] = property.default;
					}
				}
				setEditorValue(defaultValue);
			}
		}
	}, [currentSchema, userInput]);

	// Debounced submit function
	const debouncedSubmit = useCallback(
		(value: any) => {
			if (debounceTimerRef.current) {
				clearTimeout(debounceTimerRef.current);
			}
			debounceTimerRef.current = setTimeout(() => {
				if (value && userInput && queueRef.current) {
					console.log("Submitting debounced value:", value, "for", userInput);
					setDisabledBtn(true);
					queueRef.current[userInput] = value;
					socket.emit("room:worker:run");
				}
			}, 500); // 500ms debounce
		},
		[userInput],
	);

	function submitEditor() {
		if (editorValue && userInput && queueRef.current) {
			setDisabledBtn(true);
			queueRef.current[userInput] = editorValue;
			socket.emit("room:worker:run");
		}
	}

	// Handle JSONForms data changes
	const handleEditorChange = useCallback((data: any) => {
		console.log("Editor value changed:", data);
		setEditorValue(data);
	}, []);

	// Handle validation errors
	const handleValidationChange = useCallback((errors: any[]) => {
		console.log("Validation errors:", errors);
		setValidationErrors(errors);
	}, []);

	useEffect(() => {
		if (sendImmediately && editorValue && userInput && queueRef.current) {
			if (initialTrigger.current) {
				initialTrigger.current = false;
				return;
			}
			debouncedSubmit(editorValue);
		}
	}, [editorValue, sendImmediately, userInput, debouncedSubmit]);

	// Cleanup debounce timer on unmount
	useEffect(() => {
		return () => {
			if (debounceTimerRef.current) {
				clearTimeout(debounceTimerRef.current);
			}
		};
	}, []);

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
			sx={{
				position: "fixed",
				top: "50px",
				left: "50px",
				height: "100%",
				maxWidth: "35%",
				margin: 0,
				padding: 0,
				display: visible ? "block" : "none",
				overflowY: "auto",
			}}
		>
			<CardHeader
				action={
					<IconButton onClick={closeMenu} aria-label="Close">
						<Close />
					</IconButton>
				}
				title={capitalizeFirstLetter(name)}
			/>
			<CardContent sx={{ paddingBottom: 10 }}>
				<Box sx={{ display: "flex", alignItems: "center", mb: 2 }}>
					<FormControl fullWidth sx={{ mr: 2 }}>
						<InputLabel>Select Option</InputLabel>
						<Select
							value={userInput || ""}
							onChange={(e) => setUserInput(e.target.value)}
							label="Select Option"
						>
							<MenuItem value="">Select...</MenuItem>
							{Object.keys(fullSchema).map((key) => (
								<MenuItem key={key} value={key}>
									{key}
								</MenuItem>
							))}
						</Select>
					</FormControl>
					{!sendImmediately && (
						<ButtonGroup>
							<Button
								variant="contained"
								color="primary"
								onClick={submitEditor}
								startIcon={<PlayArrow />}
								disabled={disabledBtn || validationErrors.length > 0}
							>
								Submit
							</Button>
							<Button
								variant="outlined"
								color="error"
								onClick={cancelTask}
								startIcon={<Stop />}
							>
								Cancel
							</Button>
						</ButtonGroup>
					)}
				</Box>
				{validationErrors.length > 0 && (
					<Alert severity="error" sx={{ mb: 2 }}>
						<strong>Validation Errors:</strong>
						<ul>
							{validationErrors.map((error, index) => (
								<li key={index}>{error.message}</li>
							))}
						</ul>
					</Alert>
				)}
				{currentSchema && (
					<JSONFormsEditor
						schema={currentSchema}
						data={editorValue}
						onChange={handleEditorChange}
						onValidationChange={handleValidationChange}
					/>
				)}
			</CardContent>
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
			<Box
				sx={{
					position: "fixed",
					top: "50px",
					left: "0",
					height: "100%",
					width: "50px",
					display: "flex",
					flexDirection: "column",
					backgroundColor: "background.paper",
					borderRight: 1,
					borderColor: "divider",
				}}
			>
				<Tooltip title="Selection" placement="right">
					<IconButton
						color={visibleOption === "selection" ? "primary" : "default"}
						onClick={() =>
							setVisibleOption(visibleOption !== "selection" ? "selection" : "")
						}
						sx={{ margin: 1 }}
					>
						<PanTool />
					</IconButton>
				</Tooltip>
				<Tooltip title="Interaction" placement="right">
					<IconButton
						color={visibleOption === "modifier" ? "primary" : "default"}
						onClick={() =>
							setVisibleOption(visibleOption !== "modifier" ? "modifier" : "")
						}
						sx={{ margin: 1 }}
					>
						<AccountTree />
					</IconButton>
				</Tooltip>
				<Tooltip title="Settings" placement="right">
					<IconButton
						color={visibleOption === "settings" ? "primary" : "default"}
						onClick={() =>
							setVisibleOption(visibleOption !== "settings" ? "settings" : "")
						}
						sx={{ margin: 1 }}
					>
						<Settings />
					</IconButton>
				</Tooltip>
				<Tooltip title="Geometry" placement="right">
					<IconButton
						color={visibleOption === "geometry" ? "primary" : "default"}
						onClick={() =>
							setVisibleOption(visibleOption !== "geometry" ? "geometry" : "")
						}
						sx={{ margin: 1 }}
					>
						<Map />
					</IconButton>
				</Tooltip>
				<Tooltip title="Analysis" placement="right">
					<IconButton
						color={visibleOption === "analysis" ? "primary" : "default"}
						onClick={() =>
							setVisibleOption(visibleOption !== "analysis" ? "analysis" : "")
						}
						sx={{ margin: 1 }}
					>
						<BarChart />
					</IconButton>
				</Tooltip>
				<Tooltip title="View on GitHub" placement="right">
					<IconButton
						component="a"
						href="https://github.com/zincware/ZnDraw"
						target="_blank"
						sx={{ margin: 1 }}
					>
						<GitHub />
					</IconButton>
				</Tooltip>
			</Box>
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
