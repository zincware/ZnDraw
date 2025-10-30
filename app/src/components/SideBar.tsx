// components/SideBar.tsx
import { Box } from "@mui/material";
import PrimaryDrawer from "./PrimaryDrawer";
import SecondaryPanel from "./SecondaryPanel";
import GeometryPanel from "./geometry/GeometryPanel";
import SelectionsPanel from "./SelectionsPanel";
import { useFormStore } from "../formStore";
import SettingsIcon from "@mui/icons-material/Settings";
import FilterCenterFocusIcon from "@mui/icons-material/FilterCenterFocus";
import BuildIcon from "@mui/icons-material/Build";
import AnalyticsIcon from '@mui/icons-material/Analytics';
import CategoryIcon from '@mui/icons-material/Category';
import { LAYOUT_CONSTANTS, getContentHeight } from "../constants/layout";


// Hardcode the navigation items
const navItems = [
  { name: "settings", icon: <SettingsIcon />, schemaType: "settings", description: "Application settings" },
  {
    name: "selections",
    icon: <FilterCenterFocusIcon />,
    schemaType: "selections",
    description: "Selection tools and groups"
  },
  { name: "modifiers", icon: <BuildIcon />, schemaType: "modifiers", description: "Modifier tools" },
  { name: "analysis", icon: <AnalyticsIcon />, schemaType: "analysis", description: "Analysis tools" },
  { name: "geometries", icon: <CategoryIcon />, schemaType: "geometries", description: "Manage geometries" },
];

const SECONDARY_PANEL_WIDTH = 240;
const SETTINGS_PANEL_WIDTH = LAYOUT_CONSTANTS.SECONDARY_DRAWER_WIDTH;
const GEOMETRIES_PANEL_WIDTH = LAYOUT_CONSTANTS.SECONDARY_DRAWER_WIDTH;
const SELECTIONS_PANEL_WIDTH = LAYOUT_CONSTANTS.SECONDARY_DRAWER_WIDTH;
const MODIFIERS_PANEL_WIDTH = LAYOUT_CONSTANTS.SECONDARY_DRAWER_WIDTH;
const ANALYSIS_PANEL_WIDTH = LAYOUT_CONSTANTS.SECONDARY_DRAWER_WIDTH;

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

  // Determine panel width based on selected category
  const panelWidth = selectedCategory === "geometries"
    ? GEOMETRIES_PANEL_WIDTH
    : selectedCategory === "selections"
    ? SELECTIONS_PANEL_WIDTH
    : selectedCategory === "settings"
    ? SETTINGS_PANEL_WIDTH
    : selectedCategory === "modifiers"
    ? MODIFIERS_PANEL_WIDTH
    : selectedCategory === "analysis"
    ? ANALYSIS_PANEL_WIDTH
    : SECONDARY_PANEL_WIDTH;

  return (
    <>
      <PrimaryDrawer
        drawerWidth={LAYOUT_CONSTANTS.PRIMARY_DRAWER_WIDTH}
        navItems={navItems}
        selectedItem={selectedCategory}
        onItemClick={handlePanelToggle}
      />
      {selectedCategory && (
        <Box
          component="aside"
          key={selectedCategory}
          sx={{
            position: "fixed",
            left: LAYOUT_CONSTANTS.PRIMARY_DRAWER_WIDTH,
            top: LAYOUT_CONSTANTS.APPBAR_HEIGHT,
            width: panelWidth,
            height: getContentHeight(),
            zIndex: (theme) => theme.zIndex.drawer,
            borderRight: "1px solid rgba(0, 0, 0, 0.12)",
            bgcolor: "background.paper",
            boxShadow: 2,
          }}
        >
          {selectedCategory === "geometries" ? (
            <GeometryPanel />
          ) : selectedCategory === "selections" ? (
            <SelectionsPanel />
          ) : (
            <SecondaryPanel
              key={selectedCategory}
              panelTitle={selectedCategory}
            />
          )}
        </Box>
      )}
    </>
  );
};

export default SideBar;
