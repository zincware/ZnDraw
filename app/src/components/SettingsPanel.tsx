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
import {
	useSettingsSchemas,
	useSettingData,
	useUpdateSetting,
} from "../hooks/useSettings";
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

	// Fetch settings schemas
	const {
		data: schemas,
		isLoading: isLoadingSchemas,
		isError: isSchemasError,
	} = useSettingsSchemas(roomId);

	// Find selected schema
	const selectedSchema = useMemo(() => {
		if (!schemas || !selectedSettingKey) return null;
		return schemas.find((s) => s.name === selectedSettingKey);
	}, [schemas, selectedSettingKey]);

	// Fetch setting data
	const {
		data: serverData,
		isLoading: isLoadingData,
		isError: isDataError,
	} = useSettingData(roomId, selectedSettingKey || "");

	// Update setting mutation
	const { mutate: updateSettingMutation } = useUpdateSetting();

	// Sync server data to local form data
	useEffect(() => {
		if (!isLoadingData && serverData !== undefined) {
			setLocalFormData(serverData ?? {});
			ignoreFirstChangeRef.current = true;
		}
	}, [isLoadingData, serverData, selectedSettingKey]);

	// Debounced submit for auto-save
	const debouncedSubmit = useMemo(
		() =>
			debounce((data: any) => {
				if (!selectedSettingKey || !roomId || !userName) return;
				updateSettingMutation({
					roomId,
					category: selectedSettingKey,
					data: data,
				});
			}, 100),
		[selectedSettingKey, roomId, userName, updateSettingMutation],
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
		const originalSchema = selectedSchema?.schema;
		if (!originalSchema) return null;
		return injectDynamicEnums(originalSchema, undefined, geometries);
	}, [selectedSchema, geometries]);

	if (isSchemasError || isDataError) {
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
				{isLoadingSchemas ? (
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
								{schemas?.map((schema) => (
									<MenuItem key={schema.name} value={schema.name}>
										<Typography>
											{schema.name
												.replace(/_/g, " ")
												.replace(/\b\w/g, (l) => l.toUpperCase())}
										</Typography>
									</MenuItem>
								))}
							</Select>
						</FormControl>

						{selectedSettingKey && (
							<>
								{isLoadingData ? (
									<FormSkeleton />
								) : dynamicSchema ? (
									<Fade in={!isLoadingData} timeout={200}>
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
								) : null}
							</>
						)}
					</>
				)}
			</Box>
		</Box>
	);
};

export default memo(SettingsPanel);
