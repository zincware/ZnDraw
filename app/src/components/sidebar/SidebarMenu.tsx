import { Close, PlayArrow, Stop } from "@mui/icons-material";
import {
	Alert,
	Box,
	Button,
	Card,
	CardContent,
	FormControl,
	IconButton,
	InputLabel,
	MenuItem,
	Select,
	Typography,
} from "@mui/material";
import isEqual from "lodash.isequal";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as znsocket from "znsocket";
import { client, socket } from "../../socket";
import JSONFormsEditor from "../JSONFormsEditor";

const SIDEBAR_WIDTH = 50;
const SIDEBAR_TOP = 50;
const SIDEBAR_BOTTOM_MARGIN = 24; // Added for bottom spacing

interface SidebarMenuProps {
	closeMenu: () => void;
	token: string;
	name: string;
	sendImmediately: boolean;
}

const SidebarMenu = ({
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
		socket.emit("schema:refresh");
		initialTrigger.current = true;

		const con = new znsocket.Dict({ client, key: `schema:${token}:${name}` });
		const sharedCon = new znsocket.Dict({
			client,
			key: `schema:default:${name}`,
		});
		const queue = new znsocket.Dict({ client, key: `queue:${token}:${name}` });
		queueRef.current = queue;

		con.entries().then((items: any) => setSchema(Object.fromEntries(items)));
		sharedCon
			.entries()
			.then((items: any) => setSharedSchema(Object.fromEntries(items)));
		queue.length().then((length: any) => setDisabledBtn(length > 0));

		const onQueueRefresh = async () =>
			setDisabledBtn((await queue.length()) > 0);
		const onSchemaRefresh = async () =>
			setSchema(Object.fromEntries(await con.entries()));
		const onSharedSchemaRefresh = async () =>
			setSharedSchema(Object.fromEntries(await sharedCon.entries()));

		queue.onRefresh(onQueueRefresh);
		con.onRefresh(onSchemaRefresh);
		sharedCon.onRefresh(onSharedSchemaRefresh);

		return () => {
			con.offRefresh();
			sharedCon.offRefresh();
			queue.offRefresh();
		};
	}, [token, name]);

	useEffect(() => {
		if (userInput === undefined) {
			const keys = Object.keys({ ...sharedSchema, ...schema });
			if (keys.length > 0) setUserInput(keys[0]);
		}
	}, [schema, sharedSchema, userInput]);

	const fullSchema = useMemo(
		() => ({ ...sharedSchema, ...schema }),
		[sharedSchema, schema],
	);
	const currentSchema = useMemo(
		() => (userInput ? fullSchema[userInput] : null),
		[fullSchema, userInput],
	);

	useEffect(() => {
		if (currentSchema && !isEqual(currentSchema, lastSchemaRef.current)) {
			lastSchemaRef.current = currentSchema;
			if (currentSchema.properties) {
				const defaultValue: any = {};
				for (const [key, prop] of Object.entries(currentSchema.properties)) {
					if ((prop as any).default !== undefined) {
						defaultValue[key] = (prop as any).default;
					}
				}
				setEditorValue(defaultValue);
			}
		}
	}, [currentSchema, userInput]);

	const submitValue = useCallback(
		(value: any) => {
			if (value && userInput && queueRef.current) {
				setDisabledBtn(true);
				queueRef.current[userInput] = value;
				socket.emit("room:worker:run");
			}
		},
		[userInput],
	);

	const debouncedSubmit = useCallback(
		(value: any) => {
			if (debounceTimerRef.current) clearTimeout(debounceTimerRef.current);
			debounceTimerRef.current = setTimeout(() => submitValue(value), 500);
		},
		[submitValue],
	);

	const submitEditor = useCallback(
		() => submitValue(editorValue),
		[submitValue, editorValue],
	);
	const handleEditorChange = useCallback(
		(data: any) => setEditorValue(data),
		[],
	);
	const handleValidationChange = useCallback(
		(errors: any[]) => setValidationErrors(errors),
		[],
	);
	const cancelTask = useCallback(
		() => queueRef.current?.clear().then(() => setDisabledBtn(false)),
		[],
	);
	const capitalizeFirstLetter = useCallback(
		(text: string) =>
			text ? text.charAt(0).toUpperCase() + text.slice(1) : "",
		[],
	);

	useEffect(() => {
		if (sendImmediately && editorValue && userInput) {
			if (initialTrigger.current) {
				initialTrigger.current = false;
				return;
			}
			debouncedSubmit(editorValue);
		}
	}, [editorValue, sendImmediately, userInput, debouncedSubmit]);

	useEffect(
		() => () => {
			if (debounceTimerRef.current) clearTimeout(debounceTimerRef.current);
		},
		[],
	);

	return (
		<Card
			sx={{
				position: "fixed",
				top: SIDEBAR_TOP,
				left: SIDEBAR_WIDTH,
				// MODIFICATION: Calculate height to respect top and bottom margins
				height: `calc(100vh - ${SIDEBAR_TOP}px - ${SIDEBAR_BOTTOM_MARGIN}px)`,
				maxWidth: "35%",
				minWidth: "350px",
				margin: 0,
				padding: 0,
				display: "flex",
				flexDirection: "column",
				boxShadow: 3,
				borderRadius: 0,
				borderLeft: "1px solid",
				borderColor: "divider",
				backgroundColor: "background.paper",
			}}
		>
			<CardContent
				sx={{
					p: 2,
					display: "flex",
					flexDirection: "column",
					height: "100%",
					overflow: "hidden",
				}}
			>
				{/* --- START: Integrated Header and Controls Block --- */}
				<Box
					sx={{
						display: "flex",
						flexDirection: "column",
						gap: 2.5,
						mb: 2,
						p: 2,
						backgroundColor: "background.default",
						borderRadius: 2,
						border: "1px solid",
						borderColor: "divider",
					}}
				>
					{/* Header: Title and Close Button */}
					<Box
						sx={{
							display: "flex",
							justifyContent: "space-between",
							alignItems: "center",
						}}
					>
						<Typography variant="h6" component="div" sx={{ fontWeight: 600 }}>
							{capitalizeFirstLetter(name)}
						</Typography>
						<IconButton onClick={closeMenu} aria-label="Close" size="small">
							<Close />
						</IconButton>
					</Box>

					{/* Dropdown Selector */}
					<FormControl fullWidth>
						<InputLabel>Select Option</InputLabel>
						<Select
							value={userInput || ""}
							onChange={(e) => setUserInput(e.target.value)}
							label="Select Option"
							size="small"
						>
							{Object.keys(fullSchema).map((key) => (
								<MenuItem
									key={key}
									value={key}
									sx={{ textTransform: "capitalize" }}
								>
									{key}
								</MenuItem>
							))}
						</Select>
					</FormControl>

					{/* Action Buttons */}
					{!sendImmediately && (
						<Box sx={{ display: "flex", gap: 1 }}>
							<Button
								variant="contained"
								color="primary"
								onClick={submitEditor}
								startIcon={<PlayArrow />}
								disabled={disabledBtn || validationErrors.length > 0}
								size="medium"
								sx={{ flex: 1, fontWeight: 600 }}
							>
								Submit
							</Button>
							<Button
								variant="outlined"
								color="error"
								onClick={cancelTask}
								startIcon={<Stop />}
								size="medium"
								sx={{ flex: 1, fontWeight: 600 }}
							>
								Cancel
							</Button>
						</Box>
					)}
				</Box>
				{/* --- END: Integrated Header and Controls Block --- */}

				{/* --- START: Dynamic Content Area --- */}
				<Box sx={{ flex: 1, overflowY: "auto", pr: 0.5, paddingBottom: 2 }}>
					{validationErrors.length > 0 && (
						<Alert severity="error" sx={{ mb: 2 }}>
							<b>Validation Errors:</b>
							<ul style={{ margin: "8px 0 0", paddingLeft: "20px" }}>
								{validationErrors.map((error, index) => (
									<li key={index}>{error.message}</li>
								))}
							</ul>
						</Alert>
					)}

					{currentSchema && (
						<Box
							sx={{
								border: "1px solid",
								borderColor: "divider",
								borderRadius: 2,
								p: 2,
							}}
						>
							<JSONFormsEditor
								schema={currentSchema}
								data={editorValue}
								onChange={handleEditorChange}
								onValidationChange={handleValidationChange}
							/>
						</Box>
					)}
				</Box>
				{/* --- END: Dynamic Content Area --- */}
			</CardContent>
		</Card>
	);
};

export default SidebarMenu;
