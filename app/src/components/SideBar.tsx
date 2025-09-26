import { Drawer, Toolbar, Box, List, ListItem, ListItemButton, ListItemIcon, FormControl, InputLabel, Select, MenuItem, Button } from "@mui/material";
import InboxIcon from '@mui/icons-material/MoveToInbox';
import MailIcon from '@mui/icons-material/Mail';
import SaveIcon from '@mui/icons-material/Save';
import { useState } from 'react';
import Divider from '@mui/material/Divider';
import ListItemText from '@mui/material/ListItemText';
import Typography from '@mui/material/Typography';
import {
	materialCells,
	materialRenderers,
} from "@jsonforms/material-renderers";
import { JsonForms } from "@jsonforms/react";

const customRenderers = [
	...materialRenderers,
	// { tester: customRangeSliderTester, renderer: CustomRangeSlider },
	// { tester: customColorPickerTester, renderer: CustomColorPicker },
	// { tester: customSmilesRendererTester, renderer: CustomSmilesRenderer },
];

const secondaryPanelWidth = 240; // Width of the secondary panel when open

const SCHEMA = {
    'ConnectedParticles': {
        'data': {
            'properties': {},
            'title': 'ConnectedParticles',
            'type': 'object'
        },
        'ui': {}
    },
    'NoneSelection': {
        'data': {
            'properties': {},
            'title': 'NoneSelection',
            'type': 'object'
        },
        'ui': {}
    },
    'All': {
        'data': {
            'description': 'Select all atoms.',
            'properties': {},
            'title': 'All',
            'type': 'object'
        },
        'ui': {}
    },
    'Invert': {
        'data': { 'properties': {}, 'title': 'Invert', 'type': 'object' },
        'ui': {}
    },
    'Range': {
        'data': {
            'properties': {
                'start': {
                    'default': 0,
                    'description': 'Start index',
                    'title': 'Start',
                    'type': 'integer'
                },
                'end': {
                    'default': 5,
                    'description': 'End index',
                    'title': 'End',
                    'type': 'integer'
                },
                'step': {
                    'default': 1,
                    'description': 'Step size',
                    'title': 'Step',
                    'type': 'integer'
                }
            },
            'title': 'Range',
            'type': 'object'
        },
        'ui': {}
    },
    'Random': {
        'data': {
            'properties': {
                'count': {
                    'description': 'Number of atoms to select',
                    'title': 'Count',
                    'type': 'integer'
                }
            },
            'required': ['count'],
            'title': 'Random',
            'type': 'object'
        },
        'ui': {}
    },
    'IdenticalSpecies': {
        'data': {
            'properties': {},
            'title': 'IdenticalSpecies',
            'type': 'object'
        },
        'ui': {}
    },
    'Neighbour': {
        'data': {
            'description': 'Select the nth order neighbours of the selected atoms.',
            'properties': {
                'order': {
                    'default': 1,
                    'description': 'Order of neighbour',
                    'title': 'Order',
                    'type': 'integer'
                }
            },
            'title': 'Neighbour',
            'type': 'object'
        },
        'ui': {}
    }
};


const SideBar = () => {
    const drawerWidth = 60;
    const navItems = ['Inbox', 'Starred', 'Send email', 'Drafts'];

    const [openPanel, setOpenPanel] = useState('Inbox');
    const [openPanelSelection, setOpenPanelSelection] = useState<string | null>(null);
    const [formData, setFormData] = useState({});

    // 2. Handler to toggle the panel's visibility.
    const handlePanelToggle = (panelName) => {
        // If the clicked panel is already open, close it by setting state to null. Otherwise, open it.
        setOpenPanel((prev) => (prev === panelName ? null : panelName));
    };

    const handleFormChange = ({ data }) => {
        setFormData(data);
    };

    const handleSubmit = () => {
        console.log('Submitting form data:', formData);
        // Add your submit logic here
    };

    return (
        <>
            <Drawer
                variant="permanent"
                sx={{
                    width: drawerWidth,
                    flexShrink: 0,
                    [`& .MuiDrawer-paper`]: { width: drawerWidth, boxSizing: 'border-box' },
                }}
            >
                <Toolbar />
                <Box sx={{ overflow: 'auto' }}>
                    <List>
                        {navItems.map((text, index) => (
                            <ListItem key={text} disablePadding>
                                {/* The button now calls the handler passed down from MainPage */}
                                <ListItemButton
                                    onClick={() => handlePanelToggle(text)}
                                    selected={openPanel === text} // Highlights the button if its panel is open
                                    sx={{ justifyContent: 'center', p: 2 }} // Styles to center the icon
                                >
                                    <ListItemIcon sx={{ minWidth: 0 }}>
                                        {index % 2 === 0 ? <InboxIcon /> : <MailIcon />}
                                    </ListItemIcon>
                                </ListItemButton>
                            </ListItem>
                        ))}
                    </List>
                </Box>
            </Drawer>
            <Box
                component="div"
                sx={{
                    width: openPanel ? secondaryPanelWidth : 0, // Animate width based on state
                    flexShrink: 0,
                    overflowX: 'hidden', // Hide content when closing
                    transition: (theme) => theme.transitions.create('width', {
                        easing: theme.transitions.easing.sharp,
                        duration: theme.transitions.duration.enteringScreen,
                    }),
                    bgcolor: 'background.paper',
                    borderRight: openPanel ? '1px solid rgba(0, 0, 0, 0.12)' : 'none',
                }}
            >
                <Toolbar />
                {/* Render content only when a panel is open */}
                {openPanel && (
                    <Box sx={{ display: 'flex', flexDirection: 'column', height: 'calc(100vh - 64px)' }}>
                        <Typography variant="h6" sx={{ p: 2, pb: 1 }}>{openPanel}</Typography>
                        <Divider />

                        <Box sx={{ p: 2, flexGrow: 1, overflow: 'auto' }}>
                            <FormControl fullWidth sx={{ mb: 2 }}>
                                <InputLabel id="demo-simple-select-label">{openPanel}</InputLabel>
                                <Select
                                    labelId="demo-simple-select-label"
                                    id="demo-simple-select"
                                    value={openPanelSelection || ''}
                                    label={openPanel}
                                    onChange={(e) => setOpenPanelSelection(e.target.value)}
                                >
                                    {Object.keys(SCHEMA).map((item, index) => (
                                        <MenuItem key={index} value={item}>{item}</MenuItem>
                                    ))}
                                </Select>
                            </FormControl>

                            {openPanelSelection && (
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

                                    <Box sx={{ mt: 1 }}>
                                        <JsonForms
                                            schema={SCHEMA[openPanelSelection]?.data}
                                            data={formData}
                                            renderers={customRenderers}
                                            cells={materialCells}
                                            onChange={handleFormChange}
                                        />
                                    </Box>
                                </>
                            )}
                        </Box>
                    </Box>
                )}
            </Box>
        </>
    );
};

export default SideBar;