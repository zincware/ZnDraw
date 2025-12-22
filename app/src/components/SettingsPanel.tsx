/**
 * SettingsPanel - Dedicated panel for user settings management.
 *
 * Settings are per-user, per-room and use dedicated settings API endpoints
 * (not the extension system). Changes auto-save with debounce.
 */

import { useEffect, useMemo, useState, useRef, useCallback, memo } from "react";
import {
	Box,
	Typography,
	Divider,
	FormControl,
	InputLabel,
	Select,
	MenuItem,
	SelectChangeEvent,
	Fade,
} from "@mui/material";
import { JsonForms } from "@jsonforms/react";
import { materialCells } from "@jsonforms/material-renderers";
import { useFormStore } from "../formStore";
import { useSettings, useUpdateSettings } from "../hooks/useSettings";
import { useAppStore } from "../store";
import { debounce } from "lodash";
import { customRenderers, injectDynamicEnums } from "../utils/jsonforms";
import { PanelSkeleton, FormSkeleton } from "./shared/LoadingSkeletons";

const SettingsPanel = () => {
	// All hooks must be called unconditionally (React Rules of Hooks)
	const roomId = useAppStore((state) => state.roomId);
	const userName = useAppStore((state) => state.userName);
	const geometries = useAppStore((state) => state.geometries);
	const [localFormData, setLocalFormData] = useState<any>({});
	const ignoreFirstChangeRef = useRef(true);

	// Get selected setting from form store
	const { selectedExtensions, setSelectedExtension } = useFormStore();
	const selectedSettingKey = selectedExtensions["settings"] || null;

	// Fetch all settings (schema + data) in one call
	// Query is disabled when roomId is missing via enabled option
	const {
		data: settingsResponse,
		isLoading,
		isError,
	} = useSettings(roomId ?? "");

	// Update settings mutation
	const { mutate: updateSettingsMutation } = useUpdateSettings();

	// Backend always returns schema and data with defaults once loaded
	const schema = settingsResponse?.schema;
	const data = settingsResponse?.data;

	// Extract category names from schema
	const categoryNames = useMemo(() => {
		if (isLoading || !schema) return [];
		return Object.keys(schema.properties);
	}, [isLoading, schema]);

	// Get schema for selected category
	const selectedSchema = useMemo(() => {
		if (isLoading || !schema || !selectedSettingKey) return null;
		const categorySchema = schema.properties[selectedSettingKey];
		const defs = schema.$defs;
		if (!categorySchema.$ref) {
			// Include $defs for nested references
			return defs ? { ...categorySchema, $defs: defs } : categorySchema;
		}
		// Resolve $ref to get the full schema
		const refName = categorySchema.$ref.split("/").pop();
		const resolved = defs?.[refName];
		if (!resolved) return null;
		// Include $defs for nested references (e.g., CameraEnum in Camera)
		return { ...resolved, $defs: defs };
	}, [isLoading, schema, selectedSettingKey]);

	// Get data for selected category
	const serverData = useMemo(() => {
		if (isLoading || !data || !selectedSettingKey) return undefined;
		return data[selectedSettingKey];
	}, [isLoading, data, selectedSettingKey]);

	// Sync server data to local form data
	useEffect(() => {
		if (!isLoading && serverData !== undefined) {
			setLocalFormData(serverData ?? {});
			ignoreFirstChangeRef.current = true;
		}
	}, [isLoading, serverData, selectedSettingKey]);

	// Debounced submit for auto-save
	const debouncedSubmit = useMemo(
		() =>
			debounce((data: any) => {
				if (!selectedSettingKey || !roomId || !userName) return;
				updateSettingsMutation({
					roomId,
					data: { [selectedSettingKey]: data },
				});
			}, 100),
		[selectedSettingKey, roomId, userName, updateSettingsMutation],
	);

	// Cleanup debounce on unmount
	useEffect(() => {
		return () => {
			debouncedSubmit.cancel();
		};
	}, [debouncedSubmit]);

	// Handle form changes with auto-save
	const handleFormChange = useCallback(
		({ data }: { data: any }) => {
			const safeData = data ?? {};
			if (ignoreFirstChangeRef.current) {
				ignoreFirstChangeRef.current = false;
				return;
			}
			setLocalFormData(safeData);
			debouncedSubmit(safeData);
		},
		[debouncedSubmit],
	);

	// Handle setting selection change
	const handleSelectionChange = (event: SelectChangeEvent<string>) => {
		setSelectedExtension("settings", event.target.value || null);
	};

	// Create dynamic schema (for property inspector, etc.)
	const dynamicSchema = useMemo(() => {
		if (!selectedSchema) return null;
		return injectDynamicEnums(selectedSchema, undefined, geometries);
	}, [selectedSchema, geometries]);

	// Early returns AFTER all hooks are called (React Rules of Hooks)
	if (!roomId || !userName) {
		return <Typography sx={{ p: 2 }}>Joining room...</Typography>;
	}

	if (isError) {
		return (
			<Typography color="error" sx={{ p: 2 }}>
				Failed to load settings.
			</Typography>
		);
	}

	return (
		<Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
			<Typography variant="h6" sx={{ p: 2, pb: 1, flexShrink: 0 }}>
				Settings
			</Typography>
			<Divider />

			<Box sx={{ p: 2, pb: 12, flexGrow: 1, overflowY: "auto" }}>
				{isLoading ? (
					<PanelSkeleton />
				) : (
					<>
						<FormControl fullWidth sx={{ mb: 2 }}>
							<InputLabel id="settings-select-label">
								Settings Category
							</InputLabel>
							<Select
								labelId="settings-select-label"
								value={selectedSettingKey || ""}
								label="Settings Category"
								onChange={handleSelectionChange}
							>
								{categoryNames.map((name) => (
									<MenuItem key={name} value={name}>
										<Typography>
											{name
												.replace(/_/g, " ")
												.replace(/\b\w/g, (l) => l.toUpperCase())}
										</Typography>
									</MenuItem>
								))}
							</Select>
						</FormControl>

						{selectedSettingKey && dynamicSchema && (
							<Fade in={true} timeout={200}>
								<Box>
									<JsonForms
										key={selectedSettingKey}
										schema={dynamicSchema}
										data={localFormData}
										renderers={customRenderers}
										cells={materialCells}
										onChange={handleFormChange}
									/>
								</Box>
							</Fade>
						)}
					</>
				)}
			</Box>
		</Box>
	);
};

export default memo(SettingsPanel);
