import {
  Drawer,
  Box,
  List,
  ListItem,
  ListItemButton,
  ListItemIcon,
  IconButton,
} from "@mui/material";

import Tooltip from "@mui/material/Tooltip";
import GitHubIcon from "@mui/icons-material/GitHub";
import { LAYOUT_CONSTANTS, getContentHeight } from "../constants/layout";

interface NavItem {
  name: string;
  icon: React.ReactElement;
  description?: string;
  schemaType: string;
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
        [`& .MuiDrawer-paper`]: {
          width: drawerWidth,
          boxSizing: "border-box",
          height: getContentHeight(),
          top: LAYOUT_CONSTANTS.APPBAR_HEIGHT,
        },
      }}
    >
      <Box sx={{ overflow: "auto", flexGrow: 1, display: "flex", flexDirection: "column" }}>
        <List sx={{ flexGrow: 1 }}>
            {navItems.map((item) => (
            <ListItem key={item.name} disablePadding>
              <Tooltip title={item.description} placement="right">
              <ListItemButton
                onClick={() => onItemClick(item.name)}
                selected={selectedItem === item.name}
                sx={{
                  justifyContent: "center",
                  alignItems: "center",
                  p: 2
                }}
              >
                <ListItemIcon sx={{
                  minWidth: 0,
                  display: "flex",
                  justifyContent: "center",
                  color: selectedItem === item.name ? "primary.main" : "inherit"
                }}>
                  {item.icon}
                </ListItemIcon>
              </ListItemButton>
              </Tooltip>
            </ListItem>
            ))}
        </List>
      </Box>
      <Box sx={{
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        p: 2,
        borderTop: "1px solid rgba(0, 0, 0, 0.12)"
      }}>
        <Tooltip title="View on GitHub" placement="right">
          <IconButton
            component="a"
            href="https://github.com/zincware/ZnDraw"
            target="_blank"
            rel="noopener noreferrer"
            aria-label="GitHub repository"
            sx={{
              color: "text.secondary",
              p: 1
            }}
          >
            <GitHubIcon />
          </IconButton>
        </Tooltip>
      </Box>
    </Drawer>
  );
};

export default PrimaryDrawer;
