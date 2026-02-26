import AddIcon from "@mui/icons-material/Add";
import {
	Box,
	Button,
	CircularProgress,
	Divider,
	Typography,
} from "@mui/material";
import { useMemo } from "react";
import { useGeometriesList } from "../../hooks/useGeometries";
import { useAppStore } from "../../store";
import { useGeometryStore } from "../../stores/geometryStore";
import GeometryForm from "./GeometryForm";
import GeometryGrid from "./GeometryGrid";

const GeometryPanel = () => {
	// Use individual selectors to prevent unnecessary re-renders
	const roomId = useAppStore((state) => state.roomId);
	const { mode, setMode, resetForm } = useGeometryStore();

	const {
		data: geometriesData,
		isLoading,
		isError,
	} = useGeometriesList(roomId);

	const geometries = useMemo(() => {
		if (!geometriesData?.items || typeof geometriesData.items !== "object") {
			return [];
		}

		// geometriesData.items is an object with geometry keys as properties
		// Each property has {type: string, data: object} structure
		return Object.keys(geometriesData.items).map((key: string) => ({
			key,
			type: geometriesData.items[key]?.type || "unknown",
		}));
	}, [geometriesData]);

	const handleCreate = () => {
		resetForm();
		setMode("create");
	};

	if (!roomId) {
		return <Typography sx={{ p: 2 }}>Joining room...</Typography>;
	}

	if (isLoading) {
		return (
			<Box
				sx={{
					display: "flex",
					justifyContent: "center",
					alignItems: "center",
					height: "100%",
				}}
			>
				<CircularProgress />
			</Box>
		);
	}

	if (isError) {
		return (
			<Typography color="error" sx={{ p: 2 }}>
				Failed to load geometries.
			</Typography>
		);
	}

	return (
		<Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
			{mode === "list" && (
				<>
					<Typography variant="h6" sx={{ p: 2, pb: 1, flexShrink: 0 }}>
						Geometries
					</Typography>
					<Divider />
					<Box sx={{ p: 2 }}>
						<Button
							variant="contained"
							startIcon={<AddIcon />}
							onClick={handleCreate}
							fullWidth
						>
							Add Geometry
						</Button>
					</Box>
					<Box sx={{ flexGrow: 1, minHeight: 0 }}>
						<GeometryGrid geometries={geometries} />
					</Box>
				</>
			)}
			{(mode === "create" || mode === "edit") && <GeometryForm />}
		</Box>
	);
};

export default GeometryPanel;
