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
import { useSchemas, useSchemaData, useSubmitAction } from '../hooks/useSchemas';
import { useAppStore } from '../store';

interface SecondaryPanelProps {
    panelTitle: string;
}

const SecondaryPanel = ({ panelTitle }: SecondaryPanelProps) => {
    const { roomId, userId } = useAppStore();
    const [localFormData, setLocalFormData] = useState<any>({});
    // const userInteractionRef = useRef(false);
    const ignoreFirstChangeRef = useRef(true);

    if (!roomId || !userId) {
        return <Typography sx={{ p: 2 }}>Joining room...</Typography>;
    }

    const { selectedMethods, setSelectedMethod } = useFormStore();
    const selectedMethod = selectedMethods[panelTitle] || null;

    const { data: schemas, isLoading: isLoadingSchemas, isError: isSchemasError } =
        useSchemas(roomId, panelTitle);

    const {
        data: serverData,
        isLoading: isLoadingData,
        isError: isDataError
    } = useSchemaData(roomId, userId, panelTitle, selectedMethod || '');

    useEffect(() => {
        if (!isLoadingData && serverData !== undefined) {
            setLocalFormData(serverData ?? {});
            // userInteractionRef.current = false;
            ignoreFirstChangeRef.current = true; // <-- added
        }
    }, [isLoadingData, serverData, selectedMethod]);

    const handleFormChange = useCallback(({ data }: { data: any }) => {
        const safeData = data ?? {};
        if (ignoreFirstChangeRef.current) {
            ignoreFirstChangeRef.current = false;
            return; // ignore JsonForms init overwrite
        }
        setLocalFormData(safeData);
        // userInteractionRef.current = true;
    }, []);

    const { mutate: submit, isPending: isSubmitting } = useSubmitAction();

    const handleSelectionChange = (event: SelectChangeEvent<string>) => {
        setSelectedMethod(panelTitle, event.target.value || null);
    };

    const handleSubmit = () => {
        if (!selectedMethod || !roomId || !userId) return;
        submit({
            roomId,
            userId,
            action: panelTitle,
            method: selectedMethod,
            data: localFormData
        });
    };

    const currentSchema = useMemo(
        () => schemas?.[selectedMethod ?? '']?.schema,
        [schemas, selectedMethod]
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

            <Box sx={{ p: 2, flexGrow: 1, overflowY: 'auto' }}>
                <FormControl fullWidth sx={{ mb: 2 }}>
                    <InputLabel id="panel-select-label">{panelTitle} Method</InputLabel>
                    <Select
                        labelId="panel-select-label"
                        value={selectedMethod || ''}
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

                {currentSchema && (
                    <>
                        <Button
                            variant="contained"
                            startIcon={<SaveIcon />}
                            onClick={handleSubmit}
                            disabled={isSubmitting || isLoadingData}
                            fullWidth
                            color="primary"
                            sx={{ mb: 2 }}
                        >
                            {isSubmitting ? 'Running...' : 'Run Action'}
                        </Button>

                        {isLoadingData ? (
                            <CircularProgress />
                        ) : (
                            <JsonForms
                                key={selectedMethod} // remount when method changes
                                schema={currentSchema}
                                data={localFormData}
                                renderers={materialRenderers}
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
