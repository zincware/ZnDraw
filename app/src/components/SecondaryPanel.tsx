// components/SecondaryPanel.tsx
import { useMemo } from 'react';
import { Box, Typography, Divider, FormControl, InputLabel, Select, MenuItem, Button, SelectChangeEvent, CircularProgress } from '@mui/material';
import SaveIcon from '@mui/icons-material/Save';
import { JsonForms } from '@jsonforms/react';
import { materialCells, materialRenderers } from "@jsonforms/material-renderers";
import { useFormStore } from '../formStore';
import { useSchemas } from '../hooks/useSchemas';
import { useAppStore } from '../store';

interface SecondaryPanelProps {
    panelTitle: string;
}

const SecondaryPanel = ({ panelTitle }: SecondaryPanelProps) => {
    // This is correct, you've used the variable name 'roomId' as requested.
    const roomId = useAppStore(state => state.roomId);

    if (!roomId) {
        return <Typography sx={{ p: 2 }}>Joining room...</Typography>;
    }

    const { data: schemas, isLoading, isError } = useSchemas(roomId, panelTitle);

    // 🚀 CORRECTED: Select each piece of state individually from the Zustand store.
    // This prevents the creation of new objects on every render and stops the infinite loop.
    const formData = useFormStore(state => state.formData);
    const uiState = useFormStore(state => state.uiState);
    const setSelectedMethod = useFormStore(state => state.setSelectedMethod);
    const updateFormData = useFormStore(state => state.updateFormData);

    const selectedMethod = uiState.selectedMethods[panelTitle] || null;
    const compositeKey = `${panelTitle}.${selectedMethod}`;
    
    const formOptions = useMemo(() => Object.keys(schemas || {}), [schemas]);
    const currentSchema = useMemo(() => schemas?.[selectedMethod ?? '']?.schema, [schemas, selectedMethod]);
    const currentFormData = formData[compositeKey] || {};

    const handleSelectionChange = (event: SelectChangeEvent<string>) => {
        setSelectedMethod(panelTitle, event.target.value || null);
    };

    const handleFormChange = ({ data }: { data: any }) => {
        if (selectedMethod) {
            updateFormData(compositeKey, data);
        }
    };

    const handleSubmit = () => {
        console.log('Submitting data for', compositeKey, 'with data:', currentFormData);
        // Here you would call your TanStack mutation, e.g., saveSettings(currentFormData)
    };
    
    if (isLoading) {
        return (
            <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' }}>
                <CircularProgress />
            </Box>
        );
    }

    if (isError) {
        return <Typography color="error" sx={{ p: 2 }}>Failed to load schemas.</Typography>;
    }

    return (
        <Box sx={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
            <Typography variant="h6" sx={{ p: 2, pb: 1, flexShrink: 0 }}>{panelTitle}</Typography>
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
                            <MenuItem key={item} value={item}>{item}</MenuItem>
                        ))}
                    </Select>
                </FormControl>
                
                {currentSchema ? (
                    <>
                        <Button
                            variant="contained"
                            startIcon={<SaveIcon />}
                            onClick={handleSubmit}
                            fullWidth
                            color="primary"
                            sx={{ mb: 2 }}
                        >
                            Run Action
                        </Button>
                        <JsonForms
                            schema={currentSchema}
                            data={currentFormData}
                            renderers={materialRenderers}
                            cells={materialCells}
                            onChange={handleFormChange}
                        />
                    </>
                ) : (
                   selectedMethod && <Typography>Select a method to see its form.</Typography>
                )}
            </Box>
        </Box>
    );
};

export default SecondaryPanel;