// components/SecondaryPanel.tsx
import { useMemo } from 'react';
import { Box, Typography, Divider, FormControl, InputLabel, Select, MenuItem, Button, SelectChangeEvent } from '@mui/material';
import SaveIcon from '@mui/icons-material/Save';
import { JsonForms } from '@jsonforms/react';
import { materialCells, materialRenderers } from "@jsonforms/material-renderers";
import { useFormStore } from '../formStore'; // Import the store

interface SecondaryPanelProps {
    panelTitle: string;
}

const SecondaryPanel = ({ panelTitle }: SecondaryPanelProps) => {
    // Select only the state and actions relevant to this component
    const formConfigs = useFormStore(state => state.formConfigs);
    const formData = useFormStore(state => state.formData);
    const selectedMethod = useFormStore(state => state.uiState.selectedMethods[panelTitle] || null);
    const setSelectedMethod = useFormStore(state => state.setSelectedMethod);
    const updateFormData = useFormStore(state => state.updateFormData);

    // Derive values from state (this logic remains the same)
    const formOptions = Object.keys(formConfigs[panelTitle] || {});
    const compositeKey = `${panelTitle}.${selectedMethod}`;
    const currentFormData = formData[compositeKey] || {};
    const selectedSchemaData = formConfigs[panelTitle]?.[selectedMethod ?? ''];
    const currentSchema = selectedSchemaData?.schema;

    // Event handlers now call store actions
    const handleSelectionChange = (event: SelectChangeEvent<string>) => {
        setSelectedMethod(panelTitle, event.target.value || null);
    };

    const handleFormChange = ({ data }: { data: any }) => {
        if (selectedMethod) {
            updateFormData(panelTitle, selectedMethod, data);
        }
    };

    const handleSubmit = () => {
        if (!selectedMethod) return;
        // The data is already in the store, but we can access it here for submission
        console.log(`Submitting for ${compositeKey}:`, currentFormData);
        alert(`Data for ${compositeKey} submitted! Check the console.`);
        // You would place your API call here, e.g., api.submit(compositeKey, currentFormData)
    };

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

                {currentSchema && (
                    <>
                        <Button
                            variant="contained"
                            startIcon={<SaveIcon />}
                            onClick={handleSubmit}
                            fullWidth
                            color="primary"
                            sx={{ mb: 2 }}
                        >
                            Save Changes
                        </Button>
                        <JsonForms
                            schema={currentSchema}
                            data={currentFormData}
                            renderers={materialRenderers}
                            cells={materialCells}
                            onChange={handleFormChange}
                        />
                    </>
                )}
            </Box>
        </Box>
    );
};

export default SecondaryPanel;