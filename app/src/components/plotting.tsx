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

// --- Type Definitions ---

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
    token: string;
    step: number;
    setStep: (step: number) => void;
    setSelectedIds: (selectedIds: Set<number>) => void;
    setSelectedFrames: (selectedFrames: IndicesState) => void;
    onRemove: (id: number) => void;
    onDuplicate: (plotName: string | null) => void;
    onUpdatePlot: (id: number, plotName: string) => void;
    bringToFront: (id: number) => void;
    zIndex: number;
}


// --------------------------------------------------------------------------
// ## Parent Component
// Manages the existence, creation, and stacking of plot windows.
// --------------------------------------------------------------------------

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

    // Effect to add a blank window via button click
    useEffect(() => {
        if (addPlotsWindow > 0) {
            const newCardId = nextId.current++;
            setCards((prev) => [...prev, { id: newCardId, plotName: null }]);
            bringToFront(newCardId);
        }
    }, [addPlotsWindow, bringToFront]);

    // Effect to automatically open windows for new plots from the server
    useEffect(() => {
        setCards(currentCards => {
            const visiblePlotNames = new Set(currentCards.map(c => c.plotName).filter(Boolean));
            const plotsToAdd = updatedPlotsList.filter(name => !visiblePlotNames.has(name));

            if (plotsToAdd.length > 0) {
                const newCards = plotsToAdd.map(plotName => {
                    const newCardId = nextId.current++;
                    bringToFront(newCardId); 
                    return { id: newCardId, plotName };
                });
                return [...currentCards, ...newCards];
            }
            return currentCards; 
        });
    }, [updatedPlotsList, bringToFront]);

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


// --------------------------------------------------------------------------
// ## Child Component
// Manages the content and interactions of a single plot window.
// --------------------------------------------------------------------------

const PlotCard = ({
    identifier,
    initialPlotName,
    token,
    step,
    setStep,
    setSelectedIds,
    setSelectedFrames,
    onRemove,
    onDuplicate,
    onUpdatePlot,
    bringToFront,
    zIndex,
}: PlotCardProps) => {
    const [availablePlots, setAvailablePlots] = useState<string[]>([]);
    const [selectedOption, setSelectedOption] = useState<string>(initialPlotName || "");
    const [rawPlotData, setRawPlotData] = useState<any[] | undefined>(undefined);
    const [plotData, setPlotData] = useState<any[] | string | undefined>(undefined);
    const [plotType, setPlotType] = useState<string>("");
    const [plotLayout, setPlotLayout] = useState<any>(undefined);
    const [isLocked, setIsLocked] = useState<boolean>(false);
    const cardRef = useRef<HTMLDivElement>(null);

    // Effect to get the list of available plots
    useEffect(() => {
        const con = new znsocket.Dict({ client, key: `room:${token}:figures` });
        const handleRefresh = async () => setAvailablePlots(await con.keys());
        con.onRefresh(handleRefresh);
        handleRefresh();
        return () => con.offRefresh(handleRefresh);
    }, [token]);
    
    // Effect to fetch plot data when selection changes
    useEffect(() => {
        if (!selectedOption) {
            setRawPlotData(undefined); // Clear data if no plot is selected
            return;
        };
        const con = new znsocket.Dict({ client, key: `room:${token}:figures` });
        con.get(selectedOption).then((data: any) => {
            if (!data) return;
            if (data._type === "plotly.graph_objs.Figure") {
                const parsed = JSON.parse(data.value);
                setPlotType("plotly");
                setRawPlotData(parsed.data);
                setPlotLayout(parsed.layout);
            } else if (data._type === "zndraw.Figure") {
                setPlotType("zndraw.Figure");
                setPlotData(data.value.base64);
                setRawPlotData(undefined); // Ensure raw data is cleared for non-plotly plots
            }
        });
    }, [selectedOption, token]);

    // Effect to process data for the step indicator
    useEffect(() => {
        if (!rawPlotData || plotType !== "plotly") {
            if (!rawPlotData) setPlotData(undefined); // Clear plot if raw data disappears
            return;
        }

        const convertToArray = (d: any) => d?.bdata ? decodeTypedArraySpec(d) : d;
        const markerList: [number, number, string][] = [];
        const updatedPlotData = rawPlotData.map((trace) => {
            const converted = { ...trace, x: convertToArray(trace.x), y: convertToArray(trace.y), customdata: convertToArray(trace.customdata) };
            converted.customdata?.forEach((cd: any, i: number) => {
                if (cd?.[0] === step) {
                    const color = converted.line?.color || converted.marker?.color || "red";
                    markerList.push([converted.x[i], converted.y[i], color]);
                }
            });
            return converted;
        });

        if (markerList.length > 0) {
            updatedPlotData.push({
                type: "scatter", mode: "markers", name: "Current Step", showlegend: false,
                x: markerList.map(m => m[0]), y: markerList.map(m => m[1]),
                marker: { color: markerList.map(m => m[2]), size: 12, symbol: "circle", line: { color: "black", width: 2 } },
            });
        }
        setPlotData(updatedPlotData);
    }, [rawPlotData, step, plotType]);

    // --- Callbacks ---
    const handleSelectChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
        const newSelected = event.target.value;
        setSelectedOption(newSelected);
        onUpdatePlot(identifier, newSelected);
    }, [identifier, onUpdatePlot]);

    const onResize = useCallback<RndResizeCallback>((e, dir, ref) => {
        setPlotLayout((prev: any) => prev ? { ...prev, width: ref.offsetWidth, height: ref.offsetHeight - 50 } : undefined);
    }, []);

    const onPlotClick = useCallback(({ points }: { points: any[] }) => {
        if (points[0]?.customdata?.[0] !== undefined) setStep(points[0].customdata[0]);
        if (points[0]?.customdata?.[1] !== undefined) setSelectedIds(new Set([points[0].customdata[1]]));
    }, [setStep, setSelectedIds]);

    const onPlotSelected = useCallback((event: any) => {
        if (!event?.points?.length) return;
        const frames = new Set<number>(event.points.map((p: any) => p.customdata ? p.customdata[0] : p.pointIndex));
        const ids = new Set<number>(event.points.map((p: any) => p.customdata?.[1]).filter(Boolean));
        if (ids.size > 0) setSelectedIds(ids);
        setSelectedFrames({ active: true, indices: frames });
    }, [setSelectedFrames, setSelectedIds]);

    const onPlotDeselect = useCallback(() => {
        setSelectedFrames({ active: true, indices: new Set() });
    }, [setSelectedFrames]);
    
    const memoizedPlotContent = useMemo(() => (
        <CardContent sx={{ flexGrow: 1, p: "2px", overflow: 'hidden' }}>
            {plotType === "plotly" && plotData ? (
                <Plot
                    data={plotData}
                    layout={{ ...plotLayout, dragmode: 'lasso', autosize: true }}
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
    ), [plotType, plotData, plotLayout, selectedOption, onPlotClick, onPlotSelected, onPlotDeselect]);

    return (
        <Rnd
            minHeight={200} minWidth={220}
            default={{ x: 20, y: 20, width: 450, height: 350 }}
            style={{ zIndex }}
            onResize={onResize}
            disableDragging={isLocked}
            bounds="window"
            dragHandleClassName="drag-handle"
        >
            <Card sx={{ width: "100%", height: "100%", display: "flex", flexDirection: "column" }} ref={cardRef}>
                <CardHeader
                    className="drag-handle"
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
                            <BtnTooltip text={isLocked ? "Unlock" : "Lock"}><IconButton onClick={() => setIsLocked(!isLocked)} size="small">{isLocked ? <Lock /> : <LockOpen />}</IconButton></BtnTooltip>
                            <BtnTooltip text="Duplicate"><IconButton color="secondary" size="small" onClick={() => onDuplicate(selectedOption)}><ContentCopy /></IconButton></BtnTooltip>
                            <BtnTooltip text="Close"><IconButton onClick={() => onRemove(identifier)} size="small"><Close /></IconButton></BtnTooltip>
                        </Box>
                    }
                />
                {memoizedPlotContent}
            </Card>
        </Rnd>
    );
};
