import { Tooltip } from "@mui/material";
import type { TooltipProps } from "@mui/material/Tooltip";

interface BtnTooltipProps {
	text: string;
	children: React.ReactElement;
	placement?: TooltipProps["placement"];
	delayShow?: number;
	delayHide?: number;
}

export const BtnTooltip: React.FC<BtnTooltipProps> = ({
	text,
	children,
	placement = "bottom",
	delayShow = 0,
	delayHide = 100,
}) => {
	return (
		<Tooltip
			title={text}
			placement={placement}
			enterDelay={delayShow}
			leaveDelay={delayHide}
		>
			{children}
		</Tooltip>
	);
};
