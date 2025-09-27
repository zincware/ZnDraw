// components/SideBar.tsx
import { useEffect, useMemo } from 'react';
import { Box } from '@mui/material';
import PrimaryDrawer from './PrimaryDrawer';
import SecondaryPanel from './SecondaryPanel';
import { useFormStore } from '../formStore'; // Assuming your store is here


const PRIMARY_DRAWER_WIDTH = 60;
const SECONDARY_PANEL_WIDTH = 240;

const SideBar = () => {
    // Select ALL state and actions needed from the store
    // const { formConfigs, selectedForm, setSelectedForm } = useFormStore(state => ({
    //     formConfigs: state.formConfigs,
    //     selectedForm: state.uiState.selectedForm,
    //     setSelectedForm: state.setSelectedForm,
    // }));

    const formConfigs = useFormStore(state => state.formConfigs);
    const selectedForm = useFormStore(state => state.uiState.selectedForm);
    const setSelectedForm = useFormStore(state => state.setSelectedForm);

    const navItems = useMemo(() => Object.keys(formConfigs), [formConfigs]);

    // The toggle logic is now handled by the store's action
    const handlePanelToggle = (panelName: string) => {
        // If clicking the same item, it will be set to null, closing the panel
        const newSelectedForm = selectedForm === panelName ? null : panelName;
        setSelectedForm(newSelectedForm);
    };

    return (
        <>
            <PrimaryDrawer
                drawerWidth={PRIMARY_DRAWER_WIDTH}
                navItems={navItems}
                selectedItem={selectedForm}
                onItemClick={handlePanelToggle}
            />
            {/* The secondary panel's visibility is now driven by the store's state */}
            {selectedForm && (
                <Box
                    component="aside"
                    sx={{
                        position: 'fixed',
                        left: PRIMARY_DRAWER_WIDTH,
                        top: 64, // Height of AppBar
                        width: SECONDARY_PANEL_WIDTH,
                        height: 'calc(100vh - 64px)',
                        zIndex: (theme) => theme.zIndex.drawer,
                        borderRight: '1px solid rgba(0, 0, 0, 0.12)',
                        bgcolor: 'background.paper',
                        boxShadow: 2,
                    }}
                >
                    {/* The SecondaryPanel no longer needs schema data passed as props */}
                    {/* It will get everything it needs from the store */}
                    <SecondaryPanel key={selectedForm} panelTitle={selectedForm} />
                </Box>
            )}
        </>
    );
};

export default SideBar;