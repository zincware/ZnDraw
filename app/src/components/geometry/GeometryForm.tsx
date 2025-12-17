import { useEffect, useState, useMemo, useRef, useCallback } from "react";
import { debounce } from "lodash";
import {
	Box,
	Button,
	TextField,
	FormControl,
	InputLabel,
	Select,
	MenuItem,
	Typography,
	CircularProgress,
	Alert,
	SelectChangeEvent,
	Chip,
	FormControlLabel,
	Switch,
} from "@mui/material";
import ArrowBackIcon from "@mui/icons-material/ArrowBack";
import { JsonForms } from "@jsonforms/react";
import { materialCells } from "@jsonforms/material-renderers";
import { useGeometryStore } from "../../stores/geometryStore";
import {
	useGeometrySchemas,
	useCreateGeometry,
} from "../../hooks/useGeometries";
import { useAppStore } from "../../store";
import { useFrameKeys } from "../../hooks/useSchemas";
import { customRenderers, injectDynamicEnums } from "../../utils/jsonforms";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";

const GeometryForm = () => {
	// Use individual selectors to prevent unnecessary re-renders
	const roomId = useAppStore((state) => state.roomId);
	const geometries = useAppStore((state) => state.geometries);
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);
	const userName = null; // we don't want to skip when saving the current userName, so undefined
	const {
		mode,
		selectedKey,
		selectedType,
		formData,
		setMode,
		setSelectedType,
		setFormData,
		resetForm,
	} = useGeometryStore();

	const [keyInput, setKeyInput] = useState("");
	const [activeState, setActiveState] = useState(true);
	const [error, setError] = useState<string | null>(null);
	const [isSaving, setIsSaving] = useState(false);
	const [isLocked, setIsLocked] = useState(false); // Track if key/type are locked
	const isInitializedRef = useRef(false);
	const isSyncingFromZustandRef = useRef(false); // Skip auto-save during external sync

	// Fetch geometry schemas
	const {
		data: schemasData,
		isLoading: isLoadingSchemas,
		isError: isSchemasError,
	} = useGeometrySchemas(roomId);

	// Fetch frame keys for dynamic enums
	const { data: metadata, isLoading: isLoadingMetadata } = useFrameKeys(
		roomId!,
	);

	// Get geometry data from Zustand store (single source of truth)
	// This ensures transform control updates are immediately visible
	// In edit mode: use selectedKey; in create mode: use keyInput (for sync after creation)
	const geometryKey = mode === "edit" ? selectedKey : keyInput.trim();
	const zustandGeometry = geometryKey ? geometries?.[geometryKey] : null;

	// Create geometry mutation
	const { mutate: createGeometry, isPending: isCreating } = useCreateGeometry();

	// Extract schemas from the response
	const schemas = useMemo(() => {
		if (!schemasData?.schemas) return {};
		return schemasData.schemas;
	}, [schemasData]);

	const schemaOptions = useMemo(() => Object.keys(schemas), [schemas]);

	// Initialize form data for edit mode using Zustand store (single source of truth)
	useEffect(() => {
		if (mode === "edit" && selectedKey && zustandGeometry) {
			isInitializedRef.current = false;
			setKeyInput(selectedKey);
			setSelectedType(zustandGeometry.type);
			setFormData(zustandGeometry.data || {});
			setActiveState(zustandGeometry.data?.active !== false); // Default to true if undefined
			setIsLocked(true); // Lock fields in edit mode

			// Mark as initialized after a short delay to avoid saving initial data
			setTimeout(() => {
				isInitializedRef.current = true;
			}, 500);
		}
		// Only depend on selectedKey change, not zustandGeometry
		// Transform field sync is handled by the effect below
		// eslint-disable-next-line react-hooks/exhaustive-deps
	}, [mode, selectedKey, setFormData, setSelectedType]);

	// Sync all geometry data fields from Zustand when they change
	// This handles updates from transform controls and external clients while the form is open
	// Zustand is the single source of truth - form always reflects Zustand state
	// Works in both edit and create modes (for transform control updates after creation)
	useEffect(() => {
		if (!isInitializedRef.current || !zustandGeometry?.data) {
			return;
		}

		// Compare Zustand data with current form data
		const zustandData = zustandGeometry.data;
		const zustandJson = JSON.stringify(zustandData);
		const formJson = JSON.stringify(formData);

		// Only update if data has changed
		if (zustandJson !== formJson) {
			// Set flag to skip auto-save (external update already saved)
			isSyncingFromZustandRef.current = true;

			// Update form data from Zustand (single source of truth)
			setFormData({ ...zustandData });
			setActiveState(zustandData?.active !== false);

			// Reset flag after debounce window
			setTimeout(() => {
				isSyncingFromZustandRef.current = false;
			}, 300);
		}
		// eslint-disable-next-line react-hooks/exhaustive-deps
	}, [zustandGeometry?.data]);

	// Check if a geometry key already exists
	const isKeyTaken = useCallback(
		(key: string): boolean => {
			if (!geometries) return false;
			return Object.keys(geometries).includes(key.trim());
		},
		[geometries],
	);

	// Apply defaults and lock fields when both key and type are present in create mode
	const lockAndCreateGeometry = useCallback(() => {
		if (
			mode === "create" &&
			keyInput.trim() &&
			selectedType &&
			geometryDefaults &&
			!isLocked
		) {
			// Check if the key is already taken
			if (isKeyTaken(keyInput)) {
				setError(
					`A geometry with the name "${keyInput.trim()}" already exists. Please choose a different name.`,
				);
				return;
			}
			const defaultData = getGeometryWithDefaults(
				{},
				selectedType,
				geometryDefaults,
			);
			setFormData(defaultData);
			setIsLocked(true);
			isInitializedRef.current = true;
		}
	}, [
		mode,
		keyInput,
		selectedType,
		geometryDefaults,
		isLocked,
		setFormData,
		isKeyTaken,
	]);

	// Handler for when key or type field is blurred
	const handleFieldBlur = useCallback(() => {
		if (keyInput.trim() && selectedType) {
			lockAndCreateGeometry();
		}
	}, [keyInput, selectedType, lockAndCreateGeometry]);

	const handleTypeChange = (event: SelectChangeEvent<string>) => {
		const type = event.target.value;
		setSelectedType(type);
		setError(null);

		// Immediately apply Pydantic defaults when type is selected
		if (mode === "create" && geometryDefaults && type) {
			const defaultData = getGeometryWithDefaults({}, type, geometryDefaults);
			setFormData(defaultData);
		}
	};

	const handleFormChange = ({ data }: { data: any }) => {
		setFormData(data ?? {});
	};

	const saveGeometry = useCallback(() => {
		if (!roomId || !keyInput.trim() || !selectedType) {
			return;
		}

		setError(null);
		setIsSaving(true);

		createGeometry(
			{
				roomId,
				key: keyInput.trim(),
				geometryType: selectedType,
				geometryData: {
					...formData,
					active: activeState,
				},
			},
			{
				onSuccess: () => {
					setIsSaving(false);
				},
				onError: (error: any) => {
					setIsSaving(false);
					const errorMsg = error.message || "Failed to save geometry";
					setError(errorMsg);
					// Snackbar is now shown automatically by withAutoLock for lock failures
				},
			},
		);
	}, [roomId, keyInput, selectedType, formData, activeState, createGeometry]);

	// Create debounced version of save
	const debouncedSave = useMemo(
		() => debounce(saveGeometry, 250),
		[saveGeometry],
	);

	// Auto-save effect
	useEffect(() => {
		// Skip if not initialized or in create mode without key/type
		if (
			!isInitializedRef.current ||
			(mode === "create" && (!keyInput.trim() || !selectedType))
		) {
			return;
		}

		// Skip if syncing from Zustand (transform controls already saved)
		if (isSyncingFromZustandRef.current) {
			return;
		}

		debouncedSave();

		// Cleanup
		return () => {
			debouncedSave.cancel();
		};
	}, [formData, keyInput, selectedType, activeState, mode, debouncedSave]);

	const handleCancel = () => {
		resetForm();
		setError(null);
		setIsLocked(false);
		setKeyInput("");
		isInitializedRef.current = false;
	};

	// Create dynamic schema with injected metadata and geometries
	const currentSchema = useMemo(() => {
		if (!selectedType || !schemas[selectedType]) return null;
		return injectDynamicEnums(schemas[selectedType], metadata, geometries);
	}, [selectedType, schemas, metadata, geometries]);

	if (isLoadingSchemas || isLoadingMetadata) {
		return (
			<Box
				sx={{
					display: "flex",
					justifyContent: "center",
					alignItems: "center",
					height: "100%",
				}}
			>
				<CircularProgress />
			</Box>
		);
	}

	if (isSchemasError) {
		return (
			<Box sx={{ p: 2 }}>
				<Typography color="error">Failed to load geometry schemas.</Typography>
				<Button onClick={handleCancel} sx={{ mt: 2 }}>
					Back to List
				</Button>
			</Box>
		);
	}

	// Handle geometry not found (deleted by another client or doesn't exist)
	// Only show this error in edit mode - in create mode the geometry may not exist yet
	if (mode === "edit" && selectedKey && !zustandGeometry) {
		return (
			<Box sx={{ p: 2 }}>
				<Typography color="error">
					Geometry "{selectedKey}" not found. It may have been deleted.
				</Typography>
				<Button onClick={handleCancel} sx={{ mt: 2 }}>
					Back to List
				</Button>
			</Box>
		);
	}

	return (
		<Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
			<Box sx={{ p: 2, borderBottom: "1px solid rgba(0, 0, 0, 0.12)" }}>
				<Box
					sx={{
						display: "flex",
						alignItems: "center",
						justifyContent: "space-between",
						mb: 2,
					}}
				>
					<Button startIcon={<ArrowBackIcon />} onClick={handleCancel}>
						Back to List
					</Button>
					{isSaving && (
						<Chip
							label="Saving..."
							size="small"
							color="primary"
							sx={{ ml: 2 }}
						/>
					)}
				</Box>
				<Typography variant="h6">
					{mode === "create" ? "Create Geometry" : "Edit Geometry"}
				</Typography>
			</Box>

			<Box sx={{ p: 2, pb: 12, flexGrow: 1, overflowY: "auto" }}>
				{error && (
					<Alert severity="error" sx={{ mb: 2 }} onClose={() => setError(null)}>
						{error}
					</Alert>
				)}

				<TextField
					fullWidth
					label="Key"
					value={keyInput}
					onChange={(e) => setKeyInput(e.target.value)}
					onBlur={handleFieldBlur}
					disabled={isLocked}
					sx={{ mb: 2 }}
					required
					error={
						!isLocked &&
						mode === "create" &&
						keyInput.trim() !== "" &&
						isKeyTaken(keyInput)
					}
					helperText={
						!isLocked && mode === "create"
							? keyInput.trim() !== "" && isKeyTaken(keyInput)
								? `Name "${keyInput.trim()}" is already taken`
								: "Enter a unique name for this geometry"
							: ""
					}
				/>

				<FormControl fullWidth sx={{ mb: 2 }} required>
					<InputLabel>Type</InputLabel>
					<Select
						value={selectedType || ""}
						label="Type"
						onChange={handleTypeChange}
						onBlur={handleFieldBlur}
						disabled={isLocked}
					>
						{schemaOptions.map((type) => (
							<MenuItem key={type} value={type}>
								{type}
							</MenuItem>
						))}
					</Select>
				</FormControl>

				<FormControlLabel
					control={
						<Switch
							checked={activeState}
							onChange={(e) => setActiveState(e.target.checked)}
						/>
					}
					label="Active (visible in scene)"
					sx={{ mb: 2 }}
				/>

				{currentSchema && (
					<>
						<Typography variant="subtitle2" sx={{ mb: 2 }}>
							Geometry Data
						</Typography>
						<JsonForms
							schema={currentSchema}
							data={formData}
							renderers={customRenderers}
							cells={materialCells}
							onChange={handleFormChange}
						/>
					</>
				)}

				<Box sx={{ mt: 3 }}>
					<Button variant="outlined" onClick={handleCancel} fullWidth>
						Close
					</Button>
				</Box>
			</Box>
		</Box>
	);
};

export default GeometryForm;
