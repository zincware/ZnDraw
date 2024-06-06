import { OverlayTrigger, Tooltip } from "react-bootstrap";
import { Placement } from "react-bootstrap/esm/types";

interface BtnTooltipProps {
  text: string;
  children: any;
  placement?: Placement;
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
    <OverlayTrigger
      placement={placement}
      delay={{ show: delayShow, hide: delayHide }}
      overlay={<Tooltip>{text}</Tooltip>}
    >
      {children}
    </OverlayTrigger>
  );
};
