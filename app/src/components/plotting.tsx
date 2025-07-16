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
import { createContext, useCallback, useContext, useEffect, useMemo, useRef, useState } from "react";
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

interface PlottingContextType {
    token: string;
    step: number;
    setStep: (step: number) => void;
    setSelectedFrames: (state: IndicesState) => void;
    setSelectedIds: (ids: Set<number>) => void;
}

// --------------------------------------------------------------------------
// ## 1. Plotting Context
// Manages shared state to avoid passing too many props down.
// --------------------------------------------------------------------------

const PlottingContext = createContext<PlottingContextType | null>(null);

const usePlottingContext = () => {
    const context = useContext(PlottingContext);
    if (!context) {
        throw new Error("usePlottingContext must be used within a PlottingProvider");
    }
    return context;
};

// --------------------------------------------------------------------------
// ## 2. Parent Component (Provider)
// Manages the existence, creation, and stacking of plot windows.
// --------------------------------------------------------------------------

export const Plotting = (props: PlottingProps) => {
    const { token, step, setStep, setSelectedFrames, setSelectedIds, addPlotsWindow, updatedPlotsList } = props;
    const [cards, setCards] = useState<{ id: number; plotName: string | null }[]>([]);
    const [zIndices, setZIndices] = useState<{ [key: number]: number }>({});
    const [highestZ, setHighestZ] = useState(100);
    const nextId = useRef(0);

    const bringToFront = useCallback((cardId: number) => {
        setHighestZ((prevZ) => {
            const newZ = prevZ + 1;
            setZIndices((prev) => ({ ...prev, [cardId]: newZ }));
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
        setCards(prev => prev.map(card => card.id === cardId ? { ...card, plotName } : card));
    }, []);

    const removeCard = useCallback((cardId: number) => {
        setCards(prev => prev.filter(card => card.id !== cardId));
    }, []);

    const duplicateCard = useCallback((plotName: string | null) => {
        const newCardId = nextId.current++;
        setCards(prev => [...prev, { id: newCardId, plotName }]);
        bringToFront(newCardId);
    }, [bringToFront]);

    const contextValue = { token, step, setStep, setSelectedFrames, setSelectedIds };

    return (
        <PlottingContext.Provider value={contextValue}>
            {cards.map(card => (
                <PlotCard
                    key={card.id}
                    identifier={card.id}
                    initialPlotName={card.plotName}
                    isFocused={zIndices[card.id] === highestZ || cards.length <= 1}
                    zIndex={zIndices[card.id] || 100}
                    onRemove={removeCard}
                    onDuplicate={duplicateCard}
                    onUpdatePlot={updateCardPlot}
                    bringToFront={bringToFront}
                />
            ))}
        </PlottingContext.Provider>
    );
};

// --------------------------------------------------------------------------
// ## 3. Child Component (Consumer)
// Manages the content and interactions of a single plot window.
// --------------------------------------------------------------------------

interface PlotCardProps {
    identifier: number;
    initialPlotName: string | null;
    isFocused: boolean;
    zIndex: number;
    onRemove: (id: number) => void;
    onDuplicate: (plotName: string | null) => void;
    onUpdatePlot: (id: number, plotName: string) => void;
    bringToFront: (id: number) => void;
}

const PlotCard = ({
    identifier,
    initialPlotName,
    isFocused,
    zIndex,
    onRemove,
    onDuplicate,
    onUpdatePlot,
    bringToFront,
}: PlotCardProps) => {
    const { token, step, setStep, setSelectedFrames, setSelectedIds } = usePlottingContext();
    const [availablePlots, setAvailablePlots] = useState<string[]>([]);
    const [selectedOption, setSelectedOption] = useState(initialPlotName || "");
    const [plotData, setPlotData] = useState<any[] | string | undefined>();
    const [plotType, setPlotType] = useState("");
    const [plotLayout, setPlotLayout] = useState<any>();
    const [isLocked, setIsLocked] = useState(false);

    useEffect(() => {
        const con = new znsocket.Dict({ client, key: `room:${token}:figures` });
        const handleRefresh = async () => setAvailablePlots(await con.keys());
        con.onRefresh(handleRefresh);
        handleRefresh();
        return () => con.offRefresh(handleRefresh);
    }, [token]);

    useEffect(() => {
        if (!selectedOption) {
            setPlotData(undefined);
            setPlotType("");
            return;
        }
        const con = new znsocket.Dict({ client, key: `room:${token}:figures` });
        con.get(selectedOption).then((data: any) => {
            if (!data) return;
            if (data._type === "plotly.graph_objs.Figure") {
                const parsed = JSON.parse(data.value);
                const processedData = processPlotData(parsed.data, step);
                setPlotType("plotly");
                setPlotData(processedData);
                setPlotLayout(parsed.layout);
            } else if (data._type === "zndraw.Figure") {
                setPlotType("zndraw.Figure");
                setPlotData(data.value.base64);
            }
        });
    }, [selectedOption, token, step]); // Re-process data when step changes

    const onPlotClick = useCallback(({ points }: { points: any[] }) => {
        if (!isFocused) return;
        if (points[0]?.customdata?.[0] !== undefined) setStep(points[0].customdata[0]);
        if (points[0]?.customdata?.[1] !== undefined) setSelectedIds(new Set([points[0].customdata[1]]));
    }, [isFocused, setStep, setSelectedIds]);

    const onPlotSelected = useCallback((event: any) => {
        if (!isFocused || !event?.points?.length) return;
        const frames = new Set<number>(event.points.map((p: any) => p.customdata?.[0] ?? p.pointIndex));
        const ids = new Set<number>(event.points.map((p: any) => p.customdata?.[1]).filter(Boolean));
        if (ids.size > 0) setSelectedIds(ids);
        setSelectedFrames({ active: true, indices: frames });
    }, [isFocused, setSelectedFrames, setSelectedIds]);

    const onPlotDeselect = useCallback(() => {
        if (!isFocused) return;
        setSelectedFrames({ active: true, indices: new Set() });
    }, [isFocused, setSelectedFrames]);

    const memoizedPlotContent = useMemo(() => (
        <CardContent sx={{ flexGrow: 1, p: "2px", overflow: "hidden" }}>
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
    ), [plotData, plotLayout, plotType, onPlotClick, onPlotSelected, onPlotDeselect]);
    
    return (
        <Rnd
            minHeight={200} minWidth={220}
            default={{ x: 20, y: 20, width: 450, height: 350 }}
            style={{ zIndex }}
            disableDragging={isLocked}
            bounds="window"
            dragHandleClassName="drag-handle"
        >
            <Card sx={{ width: "100%", height: "100%", display: "flex", flexDirection: "column" }}>
                <CardHeader
                    className="drag-handle"
                    onMouseDown={() => bringToFront(identifier)}
                    sx={{ height: 50, flexShrink: 0, pr: 1, pl: 2, cursor: isLocked ? 'default' : 'move', '& .MuiCardHeader-content': { flexGrow: 1, minWidth: 0 } }}
                    title={
                        <TextField select value={selectedOption} onChange={(e) => { setSelectedOption(e.target.value); onUpdatePlot(identifier, e.target.value); }} size="small" variant="outlined" fullWidth>
                            {!selectedOption && <MenuItem value="" disabled>Select plot</MenuItem>}
                            {availablePlots.map((plot) => <MenuItem key={plot} value={plot}>{plot}</MenuItem>)}
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

// --- Helper Function ---
function processPlotData(rawData: any[] | undefined, step: number): any[] | undefined {
    if (!rawData) return undefined;
    
    const convert = (d: any) => d?.bdata ? decodeTypedArraySpec(d) : d;
    const markers: [number, number, string][] = [];
    
    const processedTraces = rawData.map(trace => {
        const converted = { ...trace, x: convert(trace.x), y: convert(trace.y), customdata: convert(trace.customdata) };
        converted.customdata?.forEach((cd: any, i: number) => {
            if (cd?.[0] === step) {
                const color = converted.line?.color || converted.marker?.color || "red";
                markers.push([converted.x[i], converted.y[i], color]);
            }
        });
        return converted;
    });

    if (markers.length > 0) {
        processedTraces.push({
            type: "scatter", mode: "markers", name: "Current Step", showlegend: false,
            x: markers.map(m => m[0]), y: markers.map(m => m[1]),
            marker: { color: markers.map(m => m[2]), size: 12, symbol: "circle", line: { color: "black", width: 2 } },
        });
    }

    return processedTraces;
}
