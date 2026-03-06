import CloudIcon from "@mui/icons-material/Cloud";
import GroupWorkIcon from "@mui/icons-material/GroupWork";
import HourglassEmptyIcon from "@mui/icons-material/HourglassEmpty";
import PersonIcon from "@mui/icons-material/Person";
import PublicIcon from "@mui/icons-material/Public";
import { Chip, Stack } from "@mui/material";

interface JobStatusChipsProps {
	jobName: string;
	workerCount: number;
	pendingTaskCount: number;
}

/** Display status chips for a job showing scope, worker count, and pending tasks. */
export const JobStatusChips = ({
	jobName,
	workerCount,
	pendingTaskCount,
}: JobStatusChipsProps) => {
	const rawScope = jobName.split(":")[0];

	const scopeConfig = {
		"@internal": {
			icon: <CloudIcon />,
			label: "Server",
			color: "primary" as const,
		},
		"@global": {
			icon: <PublicIcon />,
			label: "Public",
			color: "success" as const,
		},
		private: {
			icon: <PersonIcon />,
			label: "Private",
			color: "default" as const,
		},
	};

	const key =
		rawScope === "@internal" || rawScope === "@global" ? rawScope : "private";
	const { icon, label, color } = scopeConfig[key];

	return (
		<Stack
			direction="row"
			spacing={1}
			sx={{ mb: 2 }}
			flexWrap="wrap"
			useFlexGap
		>
			<Chip icon={icon} label={label} color={color} size="small" />

			{key !== "@internal" && (
				<Chip
					icon={<GroupWorkIcon />}
					label={`${workerCount} worker${workerCount !== 1 ? "s" : ""}`}
					color={workerCount > 0 ? "success" : "warning"}
					size="small"
				/>
			)}

			{pendingTaskCount > 0 && (
				<Chip
					icon={<HourglassEmptyIcon />}
					label={`${pendingTaskCount} pending`}
					color="info"
					size="small"
				/>
			)}
		</Stack>
	);
};
