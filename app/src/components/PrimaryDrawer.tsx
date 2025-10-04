import {
  Drawer,
  Toolbar,
  Box,
  List,
  ListItem,
  ListItemButton,
  ListItemIcon,
} from "@mui/material";

interface NavItem {
  name: string;
  icon: React.ReactElement;
}

interface PrimaryDrawerProps {
  drawerWidth: number;
  navItems: NavItem[];
  selectedItem: string | null;
  onItemClick: (itemName: string) => void;
}

const PrimaryDrawer = ({
  drawerWidth,
  navItems,
  selectedItem,
  onItemClick,
}: PrimaryDrawerProps) => {
  return (
    <Drawer
      variant="permanent"
      sx={{
        width: drawerWidth,
        flexShrink: 0,
        [`& .MuiDrawer-paper`]: { width: drawerWidth, boxSizing: "border-box" },
      }}
    >
      <Toolbar />
      <Box sx={{ overflow: "auto" }}>
        <List>
          {navItems.map((item) => (
            <ListItem key={item.name} disablePadding>
              <ListItemButton
                onClick={() => onItemClick(item.name)}
                selected={selectedItem === item.name}
                sx={{ justifyContent: "center", p: 2 }}
              >
                <ListItemIcon sx={{ minWidth: 0 }}>{item.icon}</ListItemIcon>
              </ListItemButton>
            </ListItem>
          ))}
        </List>
      </Box>
    </Drawer>
  );
};

export default PrimaryDrawer;
