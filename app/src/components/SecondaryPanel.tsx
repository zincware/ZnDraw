// components/SecondaryPanel.tsx
import { useEffect, useMemo, useState, useRef, useCallback } from 'react';
import {
    Box,
    Typography,
    Divider,
    FormControl,
    InputLabel,
    Select,
    MenuItem,
    Button,
    SelectChangeEvent,
    CircularProgress
} from '@mui/material';
import SaveIcon from '@mui/icons-material/Save';
import { JsonForms } from '@jsonforms/react';
import { materialCells, materialRenderers } from "@jsonforms/material-renderers";
import { useFormStore } from '../formStore';
import { useSchemas, useExtensionData, useSubmitExtension } from '../hooks/useSchemas';
import { useAppStore } from '../store';
import { ExtensionStatusChips } from './ExtensionStatusChips';
import { debounce } from 'lodash';
import CustomColorPicker, { customColorPickerTester } from './jsonforms-renderers/CustomColorPicker';
import CustomRangeSlider, { customRangeSliderTester } from './jsonforms-renderers/CustomRangeSlider';

interface SecondaryPanelProps {
    panelTitle: string;
}

const customRenderers = [
	...materialRenderers,
	{ tester: customColorPickerTester, renderer: CustomColorPicker },
    { tester: customRangeSliderTester, renderer: CustomRangeSlider },
];

const SecondaryPanel = ({ panelTitle }: SecondaryPanelProps) => {
    const { roomId, userId } = useAppStore();
    const [localFormData, setLocalFormData] = useState<any>({});
    // const userInteractionRef = useRef(false);
    const ignoreFirstChangeRef = useRef(true);

    if (!roomId || !userId) {
        return <Typography sx={{ p: 2 }}>Joining room...</Typography>;
    }

    const { selectedExtensions, setSelectedExtension } = useFormStore();
    const selectedExtension = selectedExtensions[panelTitle] || null;

    const { data: schemas, isLoading: isLoadingSchemas, isError: isSchemasError } =
        useSchemas(roomId, panelTitle);

    const {
        data: serverData,
        isLoading: isLoadingData,
        isError: isDataError
    } = useExtensionData(roomId, userId, panelTitle, selectedExtension || '');

    useEffect(() => {
        if (!isLoadingData && serverData !== undefined) {
            setLocalFormData(serverData ?? {});
            // userInteractionRef.current = false;
            ignoreFirstChangeRef.current = true; // <-- added
        }
    }, [isLoadingData, serverData, selectedExtension]);

    const { mutate: submit, isPending: isSubmitting } = useSubmitExtension();

    // Debounced submit function for settings auto-save
    const debouncedSubmit = useMemo(
        () => debounce((data: any) => {
            if (!selectedExtension || !roomId || !userId) return;
            submit({
                roomId,
                userId,
                category: panelTitle,
                extension: selectedExtension,
                data: data
            });
        }, 10),
        [selectedExtension, roomId, userId, panelTitle, submit]
    );

    // Cleanup debounce on unmount
    useEffect(() => {
        return () => {
            debouncedSubmit.cancel();
        };
    }, [debouncedSubmit]);

    const handleFormChange = useCallback(({ data }: { data: any }) => {
        const safeData = data ?? {};
        if (ignoreFirstChangeRef.current) {
            ignoreFirstChangeRef.current = false;
            return; // ignore JsonForms init overwrite
        }
        setLocalFormData(safeData);

        // Auto-submit for settings category
        if (panelTitle === 'settings') {
            debouncedSubmit(safeData);
        }
        // userInteractionRef.current = true;
    }, [panelTitle, debouncedSubmit]);

    const handleSelectionChange = (event: SelectChangeEvent<string>) => {
        setSelectedExtension(panelTitle, event.target.value || null);
    };

    const handleSubmit = () => {
        if (!selectedExtension || !roomId || !userId) return;
        submit({
            roomId,
            userId,
            category: panelTitle,
            extension: selectedExtension,
            data: localFormData
        });
    };

    const currentSchema = useMemo(
        () => schemas?.[selectedExtension ?? '']?.schema,
        [schemas, selectedExtension]
    );
    const formOptions = useMemo(() => Object.keys(schemas || {}), [schemas]);

    useEffect(() => {
        console.log("localFormData changed: ", localFormData);
    }, [localFormData]);

    if (isLoadingSchemas) {
        return (
            <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' }}>
                <CircularProgress />
            </Box>
        );
    }

    if (isSchemasError || isDataError) {
        return (
            <Typography color="error" sx={{ p: 2 }}>
                Failed to load configuration.
            </Typography>
        );
    }

    return (
        <Box sx={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
            <Typography variant="h6" sx={{ p: 2, pb: 1, flexShrink: 0 }}>
                {panelTitle}
            </Typography>
            <Divider />

            <Box sx={{ p: 2, pb: 12, flexGrow: 1, overflowY: 'auto' }}>
                <FormControl fullWidth sx={{ mb: 2 }}>
                    <InputLabel id="panel-select-label">{panelTitle} Method</InputLabel>
                    <Select
                        labelId="panel-select-label"
                        value={selectedExtension || ''}
                        label={`${panelTitle} Method`}
                        onChange={handleSelectionChange}
                    >
                        {formOptions.map((item) => (
                            <MenuItem key={item} value={item}>
                                {item}
                            </MenuItem>
                        ))}
                    </Select>
                </FormControl>

                {selectedExtension && schemas && schemas[selectedExtension] && (
                    <ExtensionStatusChips metadata={schemas[selectedExtension]} />
                )}

                {currentSchema && (
                    <>
                        {/* Only show Run Extension button for non-settings categories */}
                        {panelTitle !== 'settings' && (
                            <Button
                                variant="contained"
                                startIcon={<SaveIcon />}
                                onClick={handleSubmit}
                                disabled={isSubmitting || isLoadingData}
                                fullWidth
                                color="primary"
                                sx={{ mb: 2 }}
                            >
                                {isSubmitting ? 'Running...' : 'Run Extension'}
                            </Button>
                        )}

                        {isLoadingData ? (
                            <CircularProgress />
                        ) : (
                            <JsonForms
                                key={selectedExtension} // remount when extension changes
                                schema={currentSchema}
                                data={localFormData}
                                renderers={customRenderers}
                                cells={materialCells}
                                onChange={handleFormChange}
                            />
                        )}
                    </>
                )}
            </Box>
        </Box>
    );
};

export default SecondaryPanel;
