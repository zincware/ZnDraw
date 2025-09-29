// components/SideBar.tsx
import { Box } from '@mui/material';
import PrimaryDrawer from './PrimaryDrawer';
import SecondaryPanel from './SecondaryPanel';
import { useFormStore } from '../formStore';
import SettingsIcon from '@mui/icons-material/Settings';
import FilterCenterFocusIcon from '@mui/icons-material/FilterCenterFocus';
import BuildIcon from '@mui/icons-material/Build';

// Hardcode the navigation items
const navItems = [
    { name: 'settings', icon: <SettingsIcon />, schemaType: 'settings' },
    { name: 'selections', icon: <FilterCenterFocusIcon />, schemaType: 'selections' },
    { name: 'modifiers', icon: <BuildIcon />, schemaType: 'modifiers' },
];

const PRIMARY_DRAWER_WIDTH = 60;
const SECONDARY_PANEL_WIDTH = 240;

const SideBar = () => {
    // CORRECTED: Select state directly from the top level of the store.
    // The state object no longer has a 'uiState' property.
    const selectedForm = useFormStore(state => state.selectedForm);
    const setSelectedForm = useFormStore(state => state.setSelectedForm);

    const handlePanelToggle = (panelName: string) => {
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
                    <SecondaryPanel key={selectedForm} panelTitle={selectedForm} />
                </Box>
            )}
        </>
    );
};

export default SideBar;