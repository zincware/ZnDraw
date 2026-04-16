import { useTheme } from "@mui/material/styles";

/**
 * Writes the active MUI palette into `--mui-palette-*` CSS custom properties
 * on :root. The app uses `createTheme({ colorSchemes: { dark: true } })`
 * WITHOUT `cssVariables: true` so that `useTheme().palette.*` returns real
 * color strings (required by R3F Canvas, Plotly layout, chat, and any
 * other JS-side theme consumer). This component re-renders when the theme
 * mode changes and updates the <style> block, so pure-CSS consumers
 * (dockview-mui.css) can still read MUI tokens as CSS variables.
 */
export function MuiCssVars() {
	const theme = useTheme();
	const p = theme.palette;
	const css = `:root {
		--mui-palette-background-default: ${p.background.default};
		--mui-palette-background-paper: ${p.background.paper};
		--mui-palette-text-primary: ${p.text.primary};
		--mui-palette-text-secondary: ${p.text.secondary};
		--mui-palette-divider: ${p.divider};
		--mui-palette-action-hover: ${p.action.hover};
		--mui-palette-primary-main: ${p.primary.main};
	}`;
	return <style>{css}</style>;
}
