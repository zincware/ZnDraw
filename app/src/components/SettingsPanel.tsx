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
	const roomId = useAppStore((state) => state.roomId);
	const userName = useAppStore((state) => state.userName);
	const geometries = useAppStore((state) => state.geometries);
	const [localFormData, setLocalFormData] = useState<any>({});
	const ignoreFirstChangeRef = useRef(true);

	if (!roomId || !userName) {
		return <Typography sx={{ p: 2 }}>Joining room...</Typography>;
	}

	// Get selected setting from form store
	const { selectedExtensions, setSelectedExtension } = useFormStore();
	const selectedSettingKey = selectedExtensions["settings"] || null;

	// Fetch all settings (schema + data) in one call
	const {
		data: settingsResponse,
		isLoading: isLoading,
		isError: isError,
	} = useSettings(roomId);

	// Extract category names from schema
	const categoryNames = useMemo(() => {
		if (!settingsResponse?.schema?.properties) return [];
		return Object.keys(settingsResponse.schema.properties);
	}, [settingsResponse]);

	// Get schema for selected category
	const selectedSchema = useMemo(() => {
		if (!settingsResponse?.schema?.properties || !selectedSettingKey)
			return null;
		const categorySchema =
			settingsResponse.schema.properties[selectedSettingKey];
		if (!categorySchema?.$ref) return categorySchema;
		// Resolve $ref to get the full schema
		const refName = categorySchema.$ref.split("/").pop();
		return settingsResponse.schema.$defs?.[refName] ?? null;
	}, [settingsResponse, selectedSettingKey]);

	// Get data for selected category
	const serverData = useMemo(() => {
		if (!settingsResponse?.data || !selectedSettingKey) return undefined;
		return settingsResponse.data[selectedSettingKey];
	}, [settingsResponse, selectedSettingKey]);

	// Update settings mutation
	const { mutate: updateSettingsMutation } = useUpdateSettings();

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
