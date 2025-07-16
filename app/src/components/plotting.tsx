import { Close, ContentCopy, Lock, LockOpen } from "@mui/icons-material";
import {
    Box,
    Card,
    CardContent,
    CardHeader,
    IconButton,
    MenuItem,
    TextField,
    Typography,
} from "@mui/material";
import { decodeTypedArraySpec } from "plotly.js/src/lib/array.js";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import Plot from "react-plotly.js";
import { Rnd, type RndResizeCallback } from "react-rnd";
import * as znsocket from "znsocket";
import { client } from "../socket";
import { BtnTooltip } from "./tooltips";
import type { IndicesState } from "./utils";

// --- PROPS INTERFACES ---

interface PlottingProps {
    setStep: (step: number) => void;
    step: number;
    setSelectedFrames: (selectedFrames: IndicesState) => void;
    addPlotsWindow: number;
    setSelectedIds: (selectedIds: Set<number>) => void;
    token: string;
    updatedPlotsList: string[];
}

interface PlotCardProps {
    identifier: number;
    initialPlotName: string | null;
    token:string;
    step: number;
    updatedPlotsList: string[]; // <-- Added this prop back
    setStep: (step: number) => void;
    setSelectedIds: (selectedIds: Set<number>) => void;
    setSelectedFrames: (selectedFrames: IndicesState) => void;
    onRemove: (id: number) => void;
    onDuplicate: (plotName: string | null) => void;
    onUpdatePlot: (id: number, plotName: string) => void;
    bringToFront: (id: number) => void;
    zIndex: number;
}


// --- PARENT COMPONENT ---

export const Plotting = ({
    setStep,
    step,
    setSelectedFrames,
    addPlotsWindow,
    setSelectedIds,
    token,
    updatedPlotsList,
}: PlottingProps) => {
    const [cards, setCards] = useState<{ id: number; plotName: string | null }[]>([]);
    const [zIndices, setZIndices] = useState<{ [key: number]: number }>({});
    const [highestZ, setHighestZ] = useState(100);
    const nextId = useRef(0);

    const bringToFront = useCallback((cardId: number) => {
        setHighestZ((prevZ) => {
            const newZ = prevZ + 1;
            setZIndices((prevIndices) => ({ ...prevIndices, [cardId]: newZ }));
            return newZ;
        });
    }, []);

    useEffect(() => {
        if (addPlotsWindow > 0) {
            const newCardId = nextId.current++;
            setCards((prev) => [...prev, { id: newCardId, plotName: null }]);
            bringToFront(newCardId);
        }
    }, [addPlotsWindow, bringToFront]);

    useEffect(() => {
        const visiblePlotNames = new Set(cards.map(c => c.plotName).filter(Boolean));
        const newPlots = updatedPlotsList.filter(name => !visiblePlotNames.has(name));

        if (newPlots.length > 0) {
            const newCards = newPlots.map(plotName => {
                const newCardId = nextId.current++;
                bringToFront(newCardId);
                return { id: newCardId, plotName };
            });
            setCards((prev) => [...prev, ...newCards]);
        }
    }, [updatedPlotsList, cards, bringToFront]);

    const updateCardPlot = useCallback((cardId: number, plotName: string) => {
        setCards((prev) =>
            prev.map((card) => (card.id === cardId ? { ...card, plotName } : card)),
        );
    }, []);

    const removeCard = useCallback((cardId: number) => {
        setCards((prev) => prev.filter((card) => card.id !== cardId));
        setZIndices((prev) => {
            const newZ = { ...prev };
            delete newZ[cardId];
            return newZ;
        });
    }, []);

    const duplicateCard = useCallback((plotName: string | null) => {
        const newCardId = nextId.current++;
        setCards((prev) => [...prev, { id: newCardId, plotName }]);
        bringToFront(newCardId);
    }, [bringToFront]);

    return (
        <>
            {cards.map((card) => (
                <PlotCard
                    key={card.id}
                    identifier={card.id}
                    initialPlotName={card.plotName}
                    token={token}
                    step={step}
                    updatedPlotsList={updatedPlotsList} // <-- Pass prop
                    setStep={setStep}
                    setSelectedIds={setSelectedIds}
                    setSelectedFrames={setSelectedFrames}
                    onRemove={removeCard}
                    onDuplicate={duplicateCard}
                    onUpdatePlot={updateCardPlot}
                    bringToFront={bringToFront}
                    zIndex={zIndices[card.id] || 100}
                />
            ))}
        </>
    );
};

// No changes needed for Plotting, only for PlotCard.
// Full component provided for clarity.

const PlotCard = ({
    identifier,
    initialPlotName,
    token,
    step,
    updatedPlotsList,
    setStep,
    setSelectedIds,
    setSelectedFrames,
    onRemove,
    onDuplicate,
    onUpdatePlot,
    bringToFront,
    zIndex,
}: PlotCardProps) => {
    const [conInterface, setConInterface] = useState<znsocket.Dict | undefined>(undefined);
    const [availablePlots, setAvailablePlots] = useState<string[]>([]);
    const [selectedOption, setSelectedOption] = useState<string>(initialPlotName || "");
    const [rawPlotData, setRawPlotData] = useState<any[] | undefined>(undefined);
    const [plotData, setPlotData] = useState<any[] | string | undefined>(undefined);
    const [plotType, setPlotType] = useState<string>("");
    const [plotLayout, setPlotLayout] = useState<any>(undefined);
    const [isLocked, setIsLocked] = useState<boolean>(false);
    const cardRef = useRef<HTMLDivElement>(null);

    // Establish connection interface
    useEffect(() => {
        const con = new znsocket.Dict({ client, key: `room:${token}:figures` });
        const handleRefresh = async () => {
            const keys = await con.keys();
            setAvailablePlots(keys);
        };
        con.onRefresh(handleRefresh);
        handleRefresh();
        setConInterface(con);
        return () => con.offRefresh(handleRefresh);
    }, [token]);
    
    // Fetch data when the selected plot changes
    useEffect(() => {
        if (!conInterface || !selectedOption) return;
        
        conInterface.get(selectedOption).then((data: any) => {
            if (data === null) return;
            if (data._type === "plotly.graph_objs.Figure") {
                setPlotType("plotly");
                const parsed = JSON.parse(data.value);
                setRawPlotData(parsed.data);
                setPlotLayout(parsed.layout);
            } else if (data._type === "zndraw.Figure") {
                setPlotType("zndraw.Figure");
                setPlotData(data.value.base64);
            }
        });
    }, [conInterface, selectedOption]);
    
    // Refresh an open plot when its data changes on the server
    useEffect(() => {
        if (conInterface && updatedPlotsList.includes(selectedOption)) {
            conInterface.get(selectedOption).then((data: any) => {
                if (data === null) return;
                const parsed = JSON.parse(data.value);
                setRawPlotData(parsed.data);
                setPlotLayout(parsed.layout);
            });
        }
    }, [updatedPlotsList, conInterface, selectedOption]);

    // Process raw data to add the step indicator
    useEffect(() => {
        if (rawPlotData && plotType === "plotly") {
            const markerList: [number, number, string][] = [];
            const convertToArray = (data: any) => data?.bdata ? decodeTypedArraySpec(data) : data;

            const updatedPlotData = rawPlotData.map((dataItem) => {
                const convertedItem = { ...dataItem, x: convertToArray(dataItem.x), y: convertToArray(dataItem.y), customdata: dataItem.customdata ? convertToArray(dataItem.customdata) : dataItem.customdata };
                if (convertedItem.customdata) {
                    convertedItem.customdata.forEach((customdata, index) => {
                        if (customdata && customdata[0] === step && convertedItem.x?.[index] !== undefined && convertedItem.y?.[index] !== undefined) {
                            let color = "red";
                            if (convertedItem.line?.color) color = convertedItem.line.color;
                            else if (convertedItem.marker?.color) color = convertedItem.marker.color;
                            markerList.push([convertedItem.x[index], convertedItem.y[index], color]);
                        }
                    });
                }
                return convertedItem;
            });

            if (markerList.length > 0) {
                updatedPlotData.push({
                    type: "scatter", mode: "markers", name: "Current Step", showlegend: false,
                    x: markerList.map((marker) => marker[0]), y: markerList.map((marker) => marker[1]),
                    marker: { color: markerList.map((marker) => marker[2]), size: 12, symbol: "circle", line: { color: "black", width: 2 } },
                });
            }
            setPlotData(updatedPlotData);
        }
    }, [rawPlotData, step, plotType]);

    // --- Callbacks ---

    const handleSelectChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
        const newSelectedOption = event.target.value;
        setSelectedOption(newSelectedOption);
        onUpdatePlot(identifier, newSelectedOption);
    }, [identifier, onUpdatePlot]);

    const onResize = useCallback<RndResizeCallback>((e, dir, ref) => {
        if (cardRef.current) {
			setPlotLayout((prev: any) => prev ? { ...prev, width: ref.offsetWidth, height: ref.offsetHeight - 50 } : undefined);
        }
    }, []);

    const onPlotClick = useCallback(({ points }: { points: any[] }) => {
		console.log("Plot clicked:", points);
        if (points[0]?.customdata?.[0] !== undefined) setStep(points[0].customdata[0]);
        if (points[0]?.customdata?.[1] !== undefined) setSelectedIds(new Set([points[0].customdata[1]]));
    }, [setStep, setSelectedIds]);

    const onPlotSelected = useCallback((event: any) => {
        if (!event || !event.points || event.points.length === 0) return;
        const selectedFrames = event.points.map((point: any) => point.customdata ? point.customdata[0] : point.pointIndex);
        const selectedIds = new Set<number>(event.points.filter((point: any) => point.customdata?.[1]).map((point: any) => point.customdata[1]));
        if (selectedIds.size > 0) setSelectedIds(selectedIds);
        setSelectedFrames({ active: true, indices: new Set(selectedFrames) });
    }, [setSelectedFrames, setSelectedIds]);

    const onPlotDeselect = useCallback(() => {
        setSelectedFrames({ active: true, indices: new Set() });
    }, [setSelectedFrames]);
    
    // useMemo is still valuable to prevent re-rendering the plot when only the header is clicked.
    const memoizedPlotContent = useMemo(() => {
        return (
            <CardContent sx={{ flexGrow: 1, p: "2px", overflow: 'hidden' }}>
                {plotType === "plotly" && plotData ? (
                    <Plot
                        data={plotData}
                        layout={{
                            ...plotLayout,
                            // Ensure lasso selection is enabled
                            dragmode: 'lasso', 
                            autosize: true
                        }}
                        useResizeHandler={true}
                        style={{ width: "100%", height: "100%" }}
                        onClick={onPlotClick}
                        onSelected={onPlotSelected}
                        onDeselect={onPlotDeselect}
                    />
                ) : plotType === "zndraw.Figure" && plotData ? (
                    <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' }}>
                        <img src={`data:image/png;base64, ${plotData}`} alt="plot" style={{ maxWidth: "100%", maxHeight: "100%", objectFit: "contain" }} />
                    </Box>
                ) : (
                    <Typography variant="h6" color="text.secondary" sx={{ m: 3, textAlign: 'center' }}>
                        {selectedOption ? "Loading..." : "No plot selected."}
                    </Typography>
                )}
            </CardContent>
        );
    }, [plotType, plotData, plotLayout, selectedOption, onPlotClick, onPlotSelected, onPlotDeselect]);

    return (
        <Rnd
            minHeight={200} minWidth={220}
            default={{ x: 20, y: 20, width: 450, height: 350 }}
            style={{ zIndex }}
            // ✅ FIX: Removed onMouseDown from the main container
            onResize={onResize}
            disableDragging={isLocked}
            bounds="window"
            dragHandleClassName="drag-handle"
        >
            <Card sx={{ width: "100%", height: "100%", display: "flex", flexDirection: "column" }} ref={cardRef}>
                <CardHeader
                    className="drag-handle"
                    // ✅ FIX: Added onMouseDown here to isolate the event to the header
                    onMouseDown={() => bringToFront(identifier)}
                    sx={{ height: 50, flexShrink: 0, pr: 1, pl: 2, cursor: isLocked ? 'default' : 'move', '& .MuiCardHeader-content': { flexGrow: 1, minWidth: 0 } }}
                    title={
                        <TextField select value={selectedOption} onChange={handleSelectChange} size="small" variant="outlined" fullWidth>
                            {!selectedOption && <MenuItem value="" disabled>Select plot</MenuItem>}
                            {availablePlots.map((plot) => (
                                <MenuItem key={plot} value={plot}>{plot}</MenuItem>
                            ))}
                        </TextField>
                    }
                    action={
                        <Box sx={{ display: 'flex', alignItems: 'center', flexShrink: 0, ml: 1 }}>
                            <BtnTooltip text={isLocked ? "Unlock" : "Lock"}>
                                <IconButton onClick={() => setIsLocked(!isLocked)} size="small">
                                    {isLocked ? <Lock /> : <LockOpen />}
                                </IconButton>
                            </BtnTooltip>
                            <BtnTooltip text="Duplicate">
                                <IconButton color="secondary" size="small" onClick={() => onDuplicate(selectedOption)}>
                                    <ContentCopy />
                                </IconButton>
                            </BtnTooltip>
                            <BtnTooltip text="Close">
                                <IconButton onClick={() => onRemove(identifier)} size="small">
                                    <Close />
                                </IconButton>
                            </BtnTooltip>
                        </Box>
                    }
                />
                {memoizedPlotContent}
            </Card>
        </Rnd>
    );
};