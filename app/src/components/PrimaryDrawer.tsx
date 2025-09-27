// components/PrimaryDrawer.tsx
// No changes needed. This component is already well-designed.
import { Drawer, Toolbar, Box, List, ListItem, ListItemButton, ListItemIcon } from "@mui/material";
import InboxIcon from '@mui/icons-material/MoveToInbox';
import MailIcon from '@mui/icons-material/Mail';

interface PrimaryDrawerProps {
    drawerWidth: number;
    navItems: string[];
    selectedItem: string | null;
    onItemClick: (itemName: string) => void;
}

const PrimaryDrawer = ({ drawerWidth, navItems, selectedItem, onItemClick }: PrimaryDrawerProps) => {
    return (
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
                            <ListItemButton
                                onClick={() => onItemClick(text)}
                                selected={selectedItem === text}
                                sx={{ justifyContent: 'center', p: 2 }}
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
    );
};

export default PrimaryDrawer;