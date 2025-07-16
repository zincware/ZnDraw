import {
    AccountTree,
    BarChart,
    GitHub,
    Map,
    PanTool,
    Settings,
} from "@mui/icons-material";
import {
    Box,
    IconButton,
    Tooltip,
} from "@mui/material";
import { useCallback, useEffect, useMemo, useState } from "react";
import { socket } from "../socket";

import SidebarMenu from "./sidebar/SidebarMenu";

// Sidebar menu configuration
const SIDEBAR_MENU_CONFIG = [
    { name: "selection", icon: PanTool, title: "Selection", sendImmediately: false },
    { name: "modifier", icon: AccountTree, title: "Interaction", sendImmediately: false },
    { name: "settings", icon: Settings, title: "Settings", sendImmediately: true },
    { name: "geometry", icon: Map, title: "Geometry", sendImmediately: false },
    { name: "analysis", icon: BarChart, title: "Analysis", sendImmediately: false },
] as const;

const SIDEBAR_WIDTH = 50;
const SIDEBAR_TOP = 50;


const SidebarButton = ({ icon: Icon, title, isActive, onClick }: { icon: any; title: string; isActive: boolean; onClick: () => void; }) => (
    <Tooltip title={title} placement="right">
        <IconButton
            color={isActive ? "primary" : "default"}
            onClick={onClick}
            sx={{ 
                m: 0.5,
                width: 40, height: 40,
                borderRadius: 1.5,
                bgcolor: isActive ? "action.selected" : "transparent",
                "&:hover": { bgcolor: "action.hover" },
            }}
        >
            <Icon />
        </IconButton>
    </Tooltip>
);

function SideBar({ token }: { token: string }) {
    const [visibleOption, setVisibleOption] = useState<string>("");
    
    const toggleOption = useCallback((option: string) => setVisibleOption(prev => prev === option ? "" : option), []);
    const closeMenu = useCallback(() => setVisibleOption(""), []);

    useEffect(() => {
        if (visibleOption) socket.emit(`${visibleOption}:schema`);
    }, [visibleOption]);

    useEffect(() => {
        const handleKeyDown = (event: KeyboardEvent) => {
            if (event.key === "Escape") closeMenu();
        };
        document.addEventListener("keydown", handleKeyDown);
        return () => document.removeEventListener("keydown", handleKeyDown);
    }, [closeMenu]);

    const activeMenuConfig = useMemo(() => 
        SIDEBAR_MENU_CONFIG.find(menu => menu.name === visibleOption),
        [visibleOption]
    );

    return (
        <>
            <Box
                sx={{
                    position: "fixed",
                    top: SIDEBAR_TOP,
                    left: 0,
                    height: "calc(100% - 50px)",
                    width: SIDEBAR_WIDTH,
                    display: "flex",
                    flexDirection: "column",
                    backgroundColor: "background.paper",
                    borderRight: 1,
                    borderColor: "divider",
                    boxShadow: 1,
                    alignItems: "center",
                    py: 1,
                    gap: 0.5
                }}
            >
                {SIDEBAR_MENU_CONFIG.map(({ name, icon, title }) => (
                    <SidebarButton
                        key={name}
                        icon={icon}
                        title={title}
                        isActive={visibleOption === name}
                        onClick={() => toggleOption(name)}
                    />
                ))}
                <Tooltip title="View on GitHub" placement="right">
                    <IconButton
                        component="a"
                        href="https://github.com/zincware/ZnDraw"
                        target="_blank"
                        sx={{ m: 0.5, width: 40, height: 40, borderRadius: 1.5, color: 'text.secondary', '&:hover': { color: 'text.primary' } }}
                    >
                        <GitHub />
                    </IconButton>
                </Tooltip>
            </Box>
            
            {activeMenuConfig && (
                <SidebarMenu
                    key={activeMenuConfig.name}
                    name={activeMenuConfig.name}
                    token={token}
                    closeMenu={closeMenu}
                    sendImmediately={activeMenuConfig.sendImmediately}
                />
            )}
        </>
    );
}

export default SideBar;
