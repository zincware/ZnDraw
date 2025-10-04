// components/SideBar.tsx
import { Box } from "@mui/material";
import PrimaryDrawer from "./PrimaryDrawer";
import SecondaryPanel from "./SecondaryPanel";
import { useFormStore } from "../formStore";
import SettingsIcon from "@mui/icons-material/Settings";
import FilterCenterFocusIcon from "@mui/icons-material/FilterCenterFocus";
import BuildIcon from "@mui/icons-material/Build";
import AnalyticsIcon from '@mui/icons-material/Analytics';


// Hardcode the navigation items
const navItems = [
  { name: "settings", icon: <SettingsIcon />, schemaType: "settings", description: "Application settings" },
  {
    name: "selections",
    icon: <FilterCenterFocusIcon />,
    schemaType: "selections",
    description: "Selection tools"
  },
  { name: "modifiers", icon: <BuildIcon />, schemaType: "modifiers", description: "Modifier tools" },
  { name: "analysis", icon: <AnalyticsIcon />, schemaType: "analysis", description: "Analysis tools" },
];

const PRIMARY_DRAWER_WIDTH = 60;
const SECONDARY_PANEL_WIDTH = 240;

const SideBar = () => {
  const selectedCategory = useFormStore((state) => state.selectedCategory);
  const setSelectedCategory = useFormStore(
    (state) => state.setSelectedCategory,
  );

  const handlePanelToggle = (panelName: string) => {
    const newSelectedCategory =
      selectedCategory === panelName ? null : panelName;
    setSelectedCategory(newSelectedCategory);
  };

  return (
    <>
      <PrimaryDrawer
        drawerWidth={PRIMARY_DRAWER_WIDTH}
        navItems={navItems}
        selectedItem={selectedCategory}
        onItemClick={handlePanelToggle}
      />
      {selectedCategory && (
        <Box
          component="aside"
          sx={{
            position: "fixed",
            left: PRIMARY_DRAWER_WIDTH,
            top: 64, // Height of AppBar
            width: SECONDARY_PANEL_WIDTH,
            height: "calc(100vh - 64px)",
            zIndex: (theme) => theme.zIndex.drawer,
            borderRight: "1px solid rgba(0, 0, 0, 0.12)",
            bgcolor: "background.paper",
            boxShadow: 2,
          }}
        >
          <SecondaryPanel
            key={selectedCategory}
            panelTitle={selectedCategory}
          />
        </Box>
      )}
    </>
  );
};

export default SideBar;
