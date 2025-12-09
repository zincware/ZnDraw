/**
 * Custom renderer for 2D vertices editing in Shape geometry.
 * Provides side-by-side table and canvas views for polygon editing.
 */

import { useState, useCallback, useMemo } from "react";
import { withJsonFormsControlProps } from "@jsonforms/react";
import { rankWith, schemaMatches, ControlProps } from "@jsonforms/core";
import {
	Box,
	FormLabel,
	IconButton,
	Tooltip,
	Table,
	TableBody,
	TableCell,
	TableContainer,
	TableHead,
	TableRow,
	Paper,
	TextField,
	Button,
	Typography,
	Slider,
	Collapse,
} from "@mui/material";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";
import DeleteIcon from "@mui/icons-material/Delete";
import AddIcon from "@mui/icons-material/Add";
import ClearIcon from "@mui/icons-material/Clear";
import CenterFocusStrongIcon from "@mui/icons-material/CenterFocusStrong";
import ZoomInIcon from "@mui/icons-material/ZoomIn";
import ZoomOutIcon from "@mui/icons-material/ZoomOut";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import { Stage, Layer, Line, Circle, Group, Text } from "react-konva";

type Vertex = [number, number];

interface VerticesTableProps {
	vertices: Vertex[];
	onChange: (vertices: Vertex[]) => void;
}

interface VerticesCanvasProps {
	vertices: Vertex[];
	onChange: (vertices: Vertex[]) => void;
}

/**
 * Table component for editing vertices with X, Y columns.
 */
const VerticesTable = ({ vertices, onChange }: VerticesTableProps) => {
	const handleValueChange = (index: number, axis: 0 | 1, value: string) => {
		const numValue = parseFloat(value);
		if (isNaN(numValue)) return;

		const newVertices = [...vertices];
		newVertices[index] = [...newVertices[index]] as Vertex;
		newVertices[index][axis] = numValue;
		onChange(newVertices);
	};

	const handleDelete = (index: number) => {
		const newVertices = vertices.filter((_, i) => i !== index);
		onChange(newVertices);
	};

	const handleAdd = () => {
		// Add a new vertex at the centroid of existing vertices, or origin if empty
		if (vertices.length === 0) {
			onChange([[0, 0]]);
		} else {
			const centroidX =
				vertices.reduce((sum, v) => sum + v[0], 0) / vertices.length;
			const centroidY =
				vertices.reduce((sum, v) => sum + v[1], 0) / vertices.length;
			onChange([...vertices, [centroidX, centroidY]]);
		}
	};

	return (
		<Box>
			<TableContainer component={Paper} sx={{ maxHeight: 300 }}>
				<Table size="small" stickyHeader>
					<TableHead>
						<TableRow>
							<TableCell sx={{ width: 40 }}>#</TableCell>
							<TableCell>X</TableCell>
							<TableCell>Y</TableCell>
							<TableCell sx={{ width: 40 }}></TableCell>
						</TableRow>
					</TableHead>
					<TableBody>
						{vertices.map((vertex, index) => (
							<TableRow key={index}>
								<TableCell>{index + 1}</TableCell>
								<TableCell>
									<TextField
										type="number"
										size="small"
										value={vertex[0]}
										onChange={(e) =>
											handleValueChange(index, 0, e.target.value)
										}
										inputProps={{ step: 0.1 }}
										sx={{ width: "100%" }}
									/>
								</TableCell>
								<TableCell>
									<TextField
										type="number"
										size="small"
										value={vertex[1]}
										onChange={(e) =>
											handleValueChange(index, 1, e.target.value)
										}
										inputProps={{ step: 0.1 }}
										sx={{ width: "100%" }}
									/>
								</TableCell>
								<TableCell>
									<IconButton
										size="small"
										onClick={() => handleDelete(index)}
										disabled={vertices.length <= 3}
									>
										<DeleteIcon fontSize="small" />
									</IconButton>
								</TableCell>
							</TableRow>
						))}
					</TableBody>
				</Table>
			</TableContainer>
			<Button
				startIcon={<AddIcon />}
				onClick={handleAdd}
				size="small"
				sx={{ mt: 1 }}
			>
				Add Vertex
			</Button>
		</Box>
	);
};

/**
 * Canvas component for visual polygon editing using React Konva.
 */
const VerticesCanvas = ({ vertices, onChange }: VerticesCanvasProps) => {
	const [scale, setScale] = useState(50); // pixels per unit
	const [offset, setOffset] = useState({ x: 200, y: 150 }); // center offset
	const [hoveredVertex, setHoveredVertex] = useState<number | null>(null);

	const canvasWidth = 400;
	const canvasHeight = 300;

	// Convert data coords to canvas coords (flip Y for standard math coords)
	const toCanvas = useCallback(
		(x: number, y: number) => ({
			x: x * scale + offset.x,
			y: -y * scale + offset.y,
		}),
		[scale, offset],
	);

	// Convert canvas coords to data coords
	const toData = useCallback(
		(canvasX: number, canvasY: number) => {
			const x = (canvasX - offset.x) / scale;
			const y = -(canvasY - offset.y) / scale;
			return { x, y };
		},
		[scale, offset],
	);

	// Handle vertex drag
	const handleVertexDrag = useCallback(
		(index: number, canvasX: number, canvasY: number) => {
			const pos = toData(canvasX, canvasY);
			const newVertices = [...vertices];
			newVertices[index] = [
				Math.round(pos.x * 1000) / 1000,
				Math.round(pos.y * 1000) / 1000,
			];
			onChange(newVertices);
		},
		[vertices, onChange, toData],
	);

	// Handle canvas click to add vertex
	const handleCanvasClick = useCallback(
		(e: any) => {
			// Only add if clicking on empty space (not on a vertex)
			if (e.target === e.target.getStage()) {
				const stage = e.target.getStage();
				const pointerPos = stage.getPointerPosition();
				if (pointerPos) {
					const pos = toData(pointerPos.x, pointerPos.y);
					const newVertex: Vertex = [
						Math.round(pos.x * 1000) / 1000,
						Math.round(pos.y * 1000) / 1000,
					];
					onChange([...vertices, newVertex]);
				}
			}
		},
		[vertices, onChange, toData],
	);

	// Handle vertex double-click to delete
	const handleVertexDoubleClick = useCallback(
		(index: number) => {
			if (vertices.length > 3) {
				const newVertices = vertices.filter((_, i) => i !== index);
				onChange(newVertices);
			}
		},
		[vertices, onChange],
	);

	// Center view on vertices
	const centerView = useCallback(() => {
		if (vertices.length === 0) {
			setOffset({ x: canvasWidth / 2, y: canvasHeight / 2 });
			return;
		}

		const minX = Math.min(...vertices.map((v) => v[0]));
		const maxX = Math.max(...vertices.map((v) => v[0]));
		const minY = Math.min(...vertices.map((v) => v[1]));
		const maxY = Math.max(...vertices.map((v) => v[1]));

		const centerX = (minX + maxX) / 2;
		const centerY = (minY + maxY) / 2;

		setOffset({
			x: canvasWidth / 2 - centerX * scale,
			y: canvasHeight / 2 + centerY * scale, // flip Y
		});
	}, [vertices, scale, canvasWidth, canvasHeight]);

	// Clear all vertices
	const clearVertices = useCallback(() => {
		onChange([]);
	}, [onChange]);

	// Polygon points for Konva Line
	const polygonPoints = useMemo(() => {
		return vertices.flatMap((v) => {
			const c = toCanvas(v[0], v[1]);
			return [c.x, c.y];
		});
	}, [vertices, toCanvas]);

	return (
		<Box>
			<Box sx={{ display: "flex", gap: 1, mb: 1, alignItems: "center" }}>
				<Tooltip title="Clear all">
					<IconButton size="small" onClick={clearVertices}>
						<ClearIcon fontSize="small" />
					</IconButton>
				</Tooltip>
				<Tooltip title="Center view">
					<IconButton size="small" onClick={centerView}>
						<CenterFocusStrongIcon fontSize="small" />
					</IconButton>
				</Tooltip>
				<Tooltip title="Zoom out">
					<IconButton
						size="small"
						onClick={() => setScale((s) => Math.max(10, s - 10))}
					>
						<ZoomOutIcon fontSize="small" />
					</IconButton>
				</Tooltip>
				<Slider
					value={scale}
					onChange={(_, v) => setScale(v as number)}
					min={10}
					max={100}
					size="small"
					sx={{ width: 80, mx: 1 }}
				/>
				<Tooltip title="Zoom in">
					<IconButton
						size="small"
						onClick={() => setScale((s) => Math.min(100, s + 10))}
					>
						<ZoomInIcon fontSize="small" />
					</IconButton>
				</Tooltip>
			</Box>

			<Paper
				elevation={2}
				sx={{
					border: "1px solid",
					borderColor: "divider",
					borderRadius: 1,
					overflow: "hidden",
				}}
			>
				<Stage
					width={canvasWidth}
					height={canvasHeight}
					onClick={handleCanvasClick}
					style={{ cursor: "crosshair" }}
				>
					<Layer>
						{/* Axis lines */}
						<Line
							points={[0, offset.y, canvasWidth, offset.y]}
							stroke="#ddd"
							strokeWidth={1}
						/>
						<Line
							points={[offset.x, 0, offset.x, canvasHeight]}
							stroke="#ddd"
							strokeWidth={1}
						/>

						{/* Origin label */}
						<Text
							x={offset.x + 5}
							y={offset.y + 5}
							text="0,0"
							fontSize={10}
							fill="#999"
						/>

						{/* Polygon fill */}
						{vertices.length >= 3 && (
							<Line
								points={polygonPoints}
								closed
								fill="rgba(25, 118, 210, 0.2)"
								stroke="#1976d2"
								strokeWidth={2}
							/>
						)}

						{/* Polygon outline when less than 3 vertices */}
						{vertices.length > 0 && vertices.length < 3 && (
							<Line points={polygonPoints} stroke="#1976d2" strokeWidth={2} />
						)}

						{/* Vertex handles */}
						{vertices.map((vertex, i) => {
							const pos = toCanvas(vertex[0], vertex[1]);
							const isHovered = hoveredVertex === i;
							return (
								<Group key={i}>
									<Circle
										x={pos.x}
										y={pos.y}
										radius={isHovered ? 10 : 8}
										fill={isHovered ? "#1565c0" : "#1976d2"}
										stroke="#fff"
										strokeWidth={2}
										draggable
										onDragMove={(e) => {
											handleVertexDrag(i, e.target.x(), e.target.y());
										}}
										onMouseEnter={() => setHoveredVertex(i)}
										onMouseLeave={() => setHoveredVertex(null)}
										onDblClick={() => handleVertexDoubleClick(i)}
										style={{ cursor: "move" }}
									/>
									<Text
										x={pos.x + 12}
										y={pos.y - 6}
										text={`${i + 1}`}
										fontSize={12}
										fill="#333"
									/>
								</Group>
							);
						})}
					</Layer>
				</Stage>
			</Paper>

			<Typography
				variant="caption"
				color="text.secondary"
				sx={{ mt: 0.5, display: "block" }}
			>
				Click to add vertex | Drag to move | Double-click to delete
			</Typography>
		</Box>
	);
};

/**
 * Main Vertices2D Renderer component.
 * Provides side-by-side table and canvas for editing 2D polygon vertices.
 */
const Vertices2DRenderer = ({
	data,
	handleChange,
	path,
	label,
	schema,
	required,
}: ControlProps) => {
	const [tableExpanded, setTableExpanded] = useState(false);

	// Normalize vertices data
	const vertices: Vertex[] = useMemo(() => {
		if (!data || !Array.isArray(data)) return [];
		return data.map((v: any) => {
			if (Array.isArray(v) && v.length >= 2) {
				return [Number(v[0]), Number(v[1])] as Vertex;
			}
			return [0, 0] as Vertex;
		});
	}, [data]);

	const handleVerticesChange = useCallback(
		(newVertices: Vertex[]) => {
			handleChange(path, newVertices);
		},
		[handleChange, path],
	);

	return (
		<Box sx={{ mb: 2 }}>
			{/* Header */}
			<Box sx={{ display: "flex", alignItems: "center", gap: 1, mb: 1 }}>
				<FormLabel required={required}>{label}</FormLabel>
				{schema.description && (
					<Tooltip title={schema.description} placement="top">
						<IconButton size="small" sx={{ padding: 0.5 }}>
							<HelpOutlineIcon fontSize="small" />
						</IconButton>
					</Tooltip>
				)}
			</Box>

			{/* Canvas Section (on top) */}
			<VerticesCanvas vertices={vertices} onChange={handleVerticesChange} />

			{/* Collapsible Table Section (below canvas) */}
			<Paper
				variant="outlined"
				sx={{
					mt: 2,
					borderRadius: 1,
					overflow: "hidden",
					borderColor: tableExpanded ? "primary.main" : "divider",
					transition: "border-color 0.2s",
				}}
			>
				<Box
					sx={{
						display: "flex",
						alignItems: "center",
						justifyContent: "space-between",
						cursor: "pointer",
						px: 1.5,
						py: 1,
						bgcolor: "action.hover",
						"&:hover": {
							bgcolor: "action.selected",
						},
					}}
					onClick={() => setTableExpanded(!tableExpanded)}
				>
					<Typography variant="subtitle2">
						Table ({vertices.length} vertices)
					</Typography>
					<IconButton size="small">
						{tableExpanded ? (
							<ExpandLessIcon fontSize="small" />
						) : (
							<ExpandMoreIcon fontSize="small" />
						)}
					</IconButton>
				</Box>
				<Collapse in={tableExpanded}>
					<Box sx={{ p: 1 }}>
						<VerticesTable
							vertices={vertices}
							onChange={handleVerticesChange}
						/>
					</Box>
				</Collapse>
			</Paper>
		</Box>
	);
};

/**
 * Tester function for JSON Forms to determine when to use this renderer.
 * Matches schema properties with x-custom-type: "vertices-2d"
 */
export const vertices2DRendererTester = rankWith(
	10, // High rank to ensure it's chosen over default renderers
	schemaMatches((schema) => {
		return (schema as any)?.["x-custom-type"] === "vertices-2d";
	}),
);

export default withJsonFormsControlProps(Vertices2DRenderer);
