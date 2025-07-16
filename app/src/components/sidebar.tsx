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
	Tooltip,
} from "@mui/material";
import isEqual from "lodash.isequal";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as znsocket from "znsocket";
import { client, socket } from "../socket";
import JSONFormsEditor from "./JSONFormsEditor";

// Sidebar menu configuration
const SIDEBAR_MENU_CONFIG = [
	{ name: "selection", icon: PanTool, title: "Selection", sendImmediately: false },
	{ name: "modifier", icon: AccountTree, title: "Interaction", sendImmediately: false },
	{ name: "settings", icon: Settings, title: "Settings", sendImmediately: true },
	{ name: "geometry", icon: Map, title: "Geometry", sendImmediately: false },
	{ name: "analysis", icon: BarChart, title: "Analysis", sendImmediately: false },
] as const;

const SIDEBAR_WIDTH = 50;
const SIDEBAR_TOP = 50;

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

	// Submit handlers
	const submitValue = useCallback((value: any) => {
		if (value && userInput && queueRef.current) {
			setDisabledBtn(true);
			queueRef.current[userInput] = value;
			socket.emit("room:worker:run");
		}
	}, [userInput]);

	const debouncedSubmit = useCallback(
		(value: any) => {
			if (debounceTimerRef.current) {
				clearTimeout(debounceTimerRef.current);
			}
			debounceTimerRef.current = setTimeout(() => submitValue(value), 500);
		},
		[submitValue],
	);

	const submitEditor = useCallback(() => submitValue(editorValue), [submitValue, editorValue]);

	// Event handlers
	const handleEditorChange = useCallback((data: any) => setEditorValue(data), []);
	const handleValidationChange = useCallback((errors: any[]) => setValidationErrors(errors), []);

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

	const cancelTask = useCallback(() => {
		queueRef.current?.clear().then(() => setDisabledBtn(false));
	}, []);

	const capitalizeFirstLetter = useCallback((text: string) => {
		return text ? text.charAt(0).toUpperCase() + text.slice(1) : "";
	}, []);

	return (
		<Card
			sx={{
				position: "fixed",
				top: SIDEBAR_TOP,
				left: SIDEBAR_WIDTH,
				height: "100%",
				maxWidth: "35%",
				minWidth: "320px",
				margin: 0,
				padding: 0,
				display: visible ? "block" : "none",
				overflowY: "auto",
				boxShadow: 3,
				borderRadius: 0,
				borderLeft: "1px solid",
				borderColor: "divider",
			}}
		>
			<CardHeader
				action={
					<IconButton 
						onClick={closeMenu} 
						aria-label="Close"
						size="small"
						sx={{
							color: "text.secondary",
							"&:hover": {
								color: "text.primary",
								backgroundColor: "action.hover"
							}
						}}
					>
						<Close />
					</IconButton>
				}
				title={capitalizeFirstLetter(name)}
				sx={{
					backgroundColor: "background.paper",
					borderBottom: "1px solid",
					borderColor: "divider",
					py: 2,
					"& .MuiCardHeader-title": {
						fontSize: "1.1rem",
						fontWeight: 600,
						color: "text.primary"
					}
				}}
			/>
			<CardContent sx={{ p: 3, paddingBottom: 10 }}>
				{/* Header Section with Select and Actions */}
				<Box sx={{ 
					display: "flex", 
					flexDirection: "column", 
					gap: 2, 
					mb: 3,
					p: 2,
					backgroundColor: "background.default",
					borderRadius: 2,
					border: "1px solid",
					borderColor: "divider"
				}}>
					<FormControl fullWidth>
						<InputLabel 
							shrink 
							sx={{ 
								fontSize: "0.875rem", 
								fontWeight: 600,
								color: "text.primary" 
							}}
						>
							Select Option
						</InputLabel>
						<Select
							value={userInput || ""}
							onChange={(e) => setUserInput(e.target.value)}
							label="Select Option"
							size="small"
							sx={{
								mt: 1,
								"& .MuiSelect-select": {
									py: 1.5,
									fontSize: "0.875rem"
								}
							}}
						>
							<MenuItem value="" disabled sx={{ fontStyle: "italic", color: "text.secondary" }}>
								Choose an option...
							</MenuItem>
							{Object.keys(fullSchema).map((key) => (
								<MenuItem key={key} value={key} sx={{ textTransform: "capitalize" }}>
									{key}
								</MenuItem>
							))}
						</Select>
					</FormControl>
					
					{!sendImmediately && (
						<Box sx={{ display: "flex", gap: 1, mt: 1 }}>
							<Button
								variant="contained"
								color="primary"
								onClick={submitEditor}
								startIcon={<PlayArrow />}
								disabled={disabledBtn || validationErrors.length > 0}
								size="small"
								sx={{
									flex: 1,
									fontWeight: 600,
									textTransform: "none",
									boxShadow: 2,
									"&:hover": {
										boxShadow: 3
									}
								}}
							>
								Submit
							</Button>
							<Button
								variant="outlined"
								color="error"
								onClick={cancelTask}
								startIcon={<Stop />}
								size="small"
								sx={{
									flex: 1,
									fontWeight: 600,
									textTransform: "none",
									borderWidth: 2,
									"&:hover": {
										borderWidth: 2
									}
								}}
							>
								Cancel
							</Button>
						</Box>
					)}
				</Box>
				{validationErrors.length > 0 && (
					<Alert 
						severity="error" 
						sx={{ 
							mb: 3, 
							borderRadius: 2,
							"& .MuiAlert-message": {
								width: "100%"
							}
						}}
					>
						<Box sx={{ display: "flex", flexDirection: "column", gap: 1 }}>
							<Box sx={{ fontWeight: 600, fontSize: "0.875rem" }}>
								Validation Errors:
							</Box>
							<Box component="ul" sx={{ 
								m: 0, 
								pl: 2, 
								"& li": { 
									mb: 0.5,
									fontSize: "0.8rem"
								} 
							}}>
								{validationErrors.map((error, index) => (
									<li key={index}>{error.message}</li>
								))}
							</Box>
						</Box>
					</Alert>
				)}
				{currentSchema && (
					<Box sx={{ 
						border: "1px solid",
						borderColor: "divider",
						borderRadius: 2,
						backgroundColor: "background.paper",
						p: 2
					}}>
						<JSONFormsEditor
							schema={currentSchema}
							data={editorValue}
							onChange={handleEditorChange}
							onValidationChange={handleValidationChange}
						/>
					</Box>
				)}
			</CardContent>
		</Card>
	);
};

// Reusable sidebar button component
const SidebarButton = ({ 
	icon: Icon, 
	title, 
	isActive, 
	onClick 
}: { 
	icon: any; 
	title: string; 
	isActive: boolean; 
	onClick: () => void; 
}) => (
	<Tooltip title={title} placement="right">
		<IconButton
			color={isActive ? "primary" : "default"}
			onClick={onClick}
			sx={{ 
				margin: 0.5,
				width: 40,
				height: 40,
				borderRadius: 1,
				backgroundColor: isActive ? "action.selected" : "transparent",
				"&:hover": {
					backgroundColor: isActive ? "action.selected" : "action.hover"
				},
				transition: "all 0.2s ease-in-out"
			}}
		>
			<Icon />
		</IconButton>
	</Tooltip>
);

function SideBar({ token }: { token: string }) {
	const [visibleOption, setVisibleOption] = useState<string>("");
	
	const toggleOption = useCallback((option: string) => {
		setVisibleOption(prev => prev === option ? "" : option);
	}, []);

	const closeMenu = useCallback(() => setVisibleOption(""), []);

	// Emit schema refresh when option changes
	useEffect(() => {
		if (visibleOption) {
			socket.emit(`${visibleOption}:schema`);
		}
	}, [visibleOption]);

	// Handle escape key to close menu
	useEffect(() => {
		const handleKeyDown = (event: KeyboardEvent) => {
			if (event.key === "Escape") {
				closeMenu();
			}
		};

		document.addEventListener("keydown", handleKeyDown);
		return () => document.removeEventListener("keydown", handleKeyDown);
	}, [closeMenu]);

	return (
		<>
			<Box
				sx={{
					position: "fixed",
					top: SIDEBAR_TOP,
					left: 0,
					height: "100%",
					width: SIDEBAR_WIDTH,
					display: "flex",
					flexDirection: "column",
					backgroundColor: "background.paper",
					borderRight: 1,
					borderColor: "divider",
					boxShadow: 1,
					alignItems: "center",
					py: 1,
					gap: 0.5
				}}
			>
				{SIDEBAR_MENU_CONFIG.map(({ name, icon, title }) => (
					<SidebarButton
						key={name}
						icon={icon}
						title={title}
						isActive={visibleOption === name}
						onClick={() => toggleOption(name)}
					/>
				))}
				<Tooltip title="View on GitHub" placement="right">
					<IconButton
						component="a"
						href="https://github.com/zincware/ZnDraw"
						target="_blank"
						sx={{ 
							margin: 0.5,
							width: 40,
							height: 40,
							borderRadius: 1,
							color: "text.secondary",
							"&:hover": {
								backgroundColor: "action.hover",
								color: "text.primary"
							},
							transition: "all 0.2s ease-in-out"
						}}
					>
						<GitHub />
					</IconButton>
				</Tooltip>
			</Box>
			
			{SIDEBAR_MENU_CONFIG.map(({ name, sendImmediately }) => (
				<SidebarMenu
					key={name}
					name={name}
					visible={visibleOption === name}
					token={token}
					closeMenu={closeMenu}
					sendImmediately={sendImmediately}
				/>
			))}
		</>
	);
}

export default SideBar;
