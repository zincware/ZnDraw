import { Close, ContentCopy, Lock, LockOpen } from "@mui/icons-material";
import {
    Box,
    Button,
    Card,
    CardContent,
    CardHeader,
    IconButton,
    TextField,
    Typography,
    MenuItem,
} from "@mui/material";
import { decodeTypedArraySpec } from "plotly.js/src/lib/array.js";
import { useEffect, useRef, useState } from "react";
import Plot from "react-plotly.js";
import { Rnd, type RndResizeCallback } from "react-rnd";
import * as znsocket from "znsocket";
import { client } from "../socket";
import { BtnTooltip } from "./tooltips";
import type { IndicesState } from "./utils";

// Assuming PlottingProps and PlotsCard2Props interfaces are correctly defined in your file.

export const Plotting = ({
    setStep,
    step,
    setSelectedFrames,
    addPlotsWindow,
    setSelectedIds,
    token,
    updatedPlotsList,
}: PlottingProps) => {
    const [displayedCards, setDisplayedCards] = useState<number[]>([]);
    const [visiblePlots, setVisiblePlots] = useState<{ [key: number]: string }>(
        {},
    );

    useEffect(() => {
        if (addPlotsWindow > 0) {
            setDisplayedCards((prevCards) => {
                const newCardIndex =
                    prevCards.length > 0 ? Math.max(...prevCards) + 1 : 0;
                return [...prevCards, newCardIndex];
            });
        }
    }, [addPlotsWindow]);

    useEffect(() => {
        updatedPlotsList.forEach((plotName) => {
            if (!Object.values(visiblePlots).includes(plotName)) {
                setDisplayedCards((prevCards) => {
                    const newCardIndex =
                        prevCards.length > 0 ? Math.max(...prevCards) + 1 : 0;
                    setVisiblePlots((prev) => ({ ...prev, [newCardIndex]: plotName }));
                    return [...prevCards, newCardIndex];
                });
            }
        });
    }, [updatedPlotsList, visiblePlots]);

    return (
        <>
            {displayedCards.map((cardIndex) => (
                <PlotsCard2
                    key={cardIndex}
                    identifier={cardIndex}
                    updatedPlotsList={updatedPlotsList}
                    token={token}
                    setVisiblePlots={setVisiblePlots}
                    setDisplayedCards={setDisplayedCards}
                    visiblePlots={visiblePlots}
                    setStep={setStep}
                    setSelectedIds={setSelectedIds}
                    setSelectedFrames={setSelectedFrames}
                    step={step}
                />
            ))}
        </>
    );
};

interface PlotsCard2Props {
    updatedPlotsList: string[];
    token: string;
    setVisiblePlots: React.Dispatch<React.SetStateAction<{ [key: number]: string }>>;
    identifier: number;
    setDisplayedCards: React.Dispatch<React.SetStateAction<number[]>>;
    visiblePlots: { [key: number]: string };
    setStep: (step: number) => void;
    setSelectedIds: (selectedIds: Set<number>) => void;
    setSelectedFrames: (selectedFrames: IndicesState) => void;
    step: number;
}


const PlotsCard2 = ({
    updatedPlotsList,
    token,
    setVisiblePlots,
    identifier,
    setDisplayedCards,
    visiblePlots,
    setStep,
    setSelectedIds,
    setSelectedFrames,
    step,
}: PlotsCard2Props) => {
    const [conInterface, setConInterface] =
        useState<znsocket.Dict | undefined>(undefined);
    const [availablePlots, setAvailablePlots] = useState<string[]>([]);
    const [selectedOption, setSelectedOption] = useState<string>("");
    const [rawPlotData, setRawPlotData] = useState<any[] | undefined>(undefined);
    const [plotData, setPlotData] = useState<any[] | string | undefined>(undefined);
    const [plotType, setPlotType] = useState<string>("");
    const [plotLayout, setPlotLayout] = useState<{ [key: string]: any } | undefined>(
        undefined,
    );
    const [isLocked, setIsLocked] = useState<boolean>(false);
    const cardRef = useRef<HTMLDivElement>(null);

    useEffect(() => {
        if (selectedOption || availablePlots.length === 0) {
            return;
        }

        const plotToSelect = availablePlots.find(
            (plotName) => !Object.values(visiblePlots).includes(plotName)
        );

        if (plotToSelect) {
            setSelectedOption(plotToSelect);
        }
    }, [availablePlots, selectedOption, visiblePlots]);

    useEffect(() => {
        if (visiblePlots[identifier] && visiblePlots[identifier] !== selectedOption) {
            setSelectedOption(visiblePlots[identifier]);
        }
    }, [identifier, visiblePlots, selectedOption]);

    useEffect(() => {
        const con = new znsocket.Dict({
            client: client,
            key: `room:${token}:figures`,
        });

        const handleRefresh = async (x: any) => {
            const keys = await con.keys();
            setAvailablePlots(keys);
        };

        con.onRefresh(handleRefresh);
        con.keys().then((keys: any) => {
            setAvailablePlots(keys);
        });
        setConInterface(con);

        return () => {
            con.offRefresh(handleRefresh);
        };
    }, [token]);

    useEffect(() => {
        if (conInterface === undefined || selectedOption === "") {
            setRawPlotData(undefined);
            setPlotLayout(undefined);
            setPlotType("");
            setPlotData(undefined);
            return;
        }

        conInterface.get(selectedOption).then((data: any) => {
            if (data === null) {
                console.warn(`No data found for plot: ${selectedOption}`);
                setRawPlotData(undefined);
                setPlotLayout(undefined);
                setPlotType("");
                setPlotData(undefined);
                return;
            }
            if (data._type === "plotly.graph_objs.Figure") {
                setPlotType("plotly");
                const parsedData = JSON.parse(data.value);
                setRawPlotData(parsedData.data);

                // --- START PLOTLY MARGIN ADJUSTMENT ---
                const initialLayout = parsedData.layout;
                setPlotLayout({
                    ...initialLayout,
                    autosize: true, // Keep autosize for container resizing
                    margin: {
                        l: 20, // left margin
                        r: 20, // right margin
                        b: 20, // bottom margin
                        t: 20, // top margin
                        pad: 0 // padding between plot area and legend/color bar
                    },
                    // You can also use automargin if you have elements that might get cut off
                    // automargin: true,
                    // If you want more control over specific axes, you can set automargin per axis
                    // xaxis: { automargin: true },
                    // yaxis: { automargin: true },
                });
                // --- END PLOTLY MARGIN ADJUSTMENT ---

            } else if (data._type === "zndraw.Figure") {
                setPlotType("zndraw.Figure");
                setPlotData(data.value.base64);
            } else {
                console.warn(`Unknown plot type: ${data._type} for plot: ${selectedOption}`);
                setRawPlotData(undefined);
                setPlotLayout(undefined);
                setPlotType("");
                setPlotData(undefined);
            }
        }).catch(error => {
            console.error(`Error fetching plot data for ${selectedOption}:`, error);
            setRawPlotData(undefined);
            setPlotLayout(undefined);
            setPlotType("");
            setPlotData(undefined);
        });
    }, [conInterface, selectedOption]);

    useEffect(() => {
        setVisiblePlots((prev) => ({ ...prev, [identifier]: selectedOption }));
    }, [selectedOption, identifier, setVisiblePlots]);

    useEffect(() => {
        if (updatedPlotsList.includes(selectedOption) && conInterface) {
            conInterface.get(selectedOption).then((data: any) => {
                if (data === null) {
                    return;
                }
                if (data._type === "plotly.graph_objs.Figure") {
                    const parsedData = JSON.parse(data.value);
                    setRawPlotData(parsedData.data);
                    
                    // Re-apply margin settings on update
                    const updatedLayout = parsedData.layout;
                    setPlotLayout({
                        ...updatedLayout,
                        autosize: true,
                        margin: {
                            l: 20,
                            r: 20,
                            b: 20,
                            t: 20,
                            pad: 0
                        },
                        // automargin: true, // Keep this consistent if you use it above
                    });
                }
            }).catch(error => {
                console.error(`Error updating plot data for ${selectedOption}:`, error);
            });
        }
    }, [updatedPlotsList, selectedOption, conInterface]);

    useEffect(() => {
        if (rawPlotData && plotType === "plotly") {
            const markerList: [number, number, string][] = [];

            const convertToArray = (data: any) => {
                if (
                    data &&
                    typeof data === "object" &&
                    "bdata" in data &&
                    "dtype" in data
                ) {
                    return decodeTypedArraySpec(data);
                }
                return data;
            };

            const updatedPlotData = rawPlotData.map((dataItem) => {
                const convertedItem = {
                    ...dataItem,
                    x: convertToArray(dataItem.x),
                    y: convertToArray(dataItem.y),
                    customdata: dataItem.customdata
                        ? convertToArray(dataItem.customdata)
                        : dataItem.customdata,
                };

                if (convertedItem.customdata) {
                    convertedItem.customdata.forEach((customdata, index) => {
                        if (
                            Array.isArray(customdata) &&
                            customdata.length > 0 &&
                            customdata[0] === step &&
                            convertedItem.x?.[index] !== undefined &&
                    convertedItem.y?.[index] !== undefined
                        ) {
                            const xPosition = convertedItem.x[index];
                            const yPosition = convertedItem.y[index];

                            let color = "red";
                            if (convertedItem.line?.color) {
                                color = convertedItem.line.color;
                            } else if (convertedItem.marker?.color) {
                                color = convertedItem.marker.color;
                            }

                            markerList.push([xPosition, yPosition, color]);
                        }
                    });
                }

                return convertedItem;
            });

            if (markerList.length > 0) {
                updatedPlotData.push({
                    type: "scatter",
                    mode: "markers",
                    name: "Current Step",
                    showlegend: false,
                    x: markerList.map((marker) => marker[0]),
                    y: markerList.map((marker) => marker[1]),
                    marker: {
                        color: markerList.map((marker) => marker[2]),
                        size: 12,
                        symbol: "circle",
                        line: {
                            color: "black",
                            width: 2,
                        },
                    },
                });
            }
            setPlotData(updatedPlotData);
        } else if (rawPlotData === undefined && plotType === "plotly") {
            setPlotData(undefined);
        }
    }, [rawPlotData, step, plotType]);

    const handleSelectChange = (event: React.ChangeEvent<HTMLInputElement>) => {
        const newSelectedOption = event.target.value;
        setSelectedOption(newSelectedOption);
    };

    const closeThisCard = () => {
        setDisplayedCards((prevCards) =>
            prevCards.filter((card) => card !== identifier),
        );
        setVisiblePlots((prevVisiblePlots) => {
            const newVisiblePlots = { ...prevVisiblePlots };
            delete newVisiblePlots[identifier];
            return newVisiblePlots;
        });
    };

    const onResize: RndResizeCallback = (e, direction, ref, delta, position) => {
        if (cardRef.current && plotLayout) { // Check plotLayout before updating
            setPlotLayout((prev) => {
                if (prev) {
                    const headerHeight = 50;
                    const contentPadding = 16;
                    // These calculations might need slight adjustments based on border/shadows if they add to the size
                    return {
                        ...prev,
                        width: ref.offsetWidth - (contentPadding * 2),
                        height: ref.offsetHeight - headerHeight - (contentPadding * 2),
                    };
                }
                return prev;
            });
        }
    };

    const onPlotClick = ({ points }: { points: any[] }) => {
        if (points && points.length > 0) {
            const firstPoint = points[0];
            if (firstPoint?.customdata?.[0] !== undefined) {
                setStep(firstPoint.customdata[0]);
            }
            if (firstPoint?.customdata?.[1] !== undefined) {
                setSelectedIds(new Set([firstPoint.customdata[1]]));
            }
        }
    };

    const addAnotherCard = () => {
        setDisplayedCards((prevCards) => {
            const newCardIndex =
                prevCards.length > 0 ? Math.max(...prevCards) + 1 : 0;
            return [...prevCards, newCardIndex];
        });
    };

    const onPlotSelected = (event: any) => {
        if (!event || !event.points || event.points.length === 0) {
            setSelectedFrames({ active: true, indices: new Set() });
            return;
        }

        const selectedFramesIndices = new Set<number>();
        const selectedObjectIds = new Set<number>();

        event.points.forEach((point: any) => {
            if (point.customdata?.[0] !== undefined) {
                selectedFramesIndices.add(point.customdata[0]);
            }
            if (point.customdata?.[1] !== undefined) {
                selectedObjectIds.add(point.customdata[1]);
            }
        });

        if (selectedObjectIds.size > 0) {
            setSelectedIds(selectedObjectIds);
        }

        setSelectedFrames({
            active: true,
            indices: selectedFramesIndices,
        });
    };

    const onPlotDeselect = () => {
        setSelectedFrames({
            active: true,
            indices: new Set(),
        });
    };

    return (
        <Rnd
            minHeight={200}
            minWidth={220}
            default={{
                x: 0,
                y: 0,
                width: 400,
                height: 300,
            }}
            onResize={onResize}
            disableDragging={isLocked}
            dragGrid={[1, 1]}
            resizeGrid={[1, 1]}
            bounds="window"
        >
            <Card
                style={{
                    width: "100%",
                    height: "100%",
                    display: "flex",
                    flexDirection: "column",
                }}
                ref={cardRef}
            >
                <CardHeader
                    sx={{
                        height: 50,
                        display: "flex",
                        alignItems: "center",
                        flexShrink: 0,
                        pr: 1,
                        pl: 2,
                        '& .MuiCardHeader-content': {
                            flexGrow: 1,
                            minWidth: 0,
                        }
                    }}
                    title={
                        <TextField
                            select
                            value={selectedOption}
                            onChange={handleSelectChange}
                            size="small"
                            variant="outlined"
                            fullWidth
                            sx={{
                                flexGrow: 1,
                            }}
                        >
                            {!selectedOption && (
                                <MenuItem value="" disabled>
                                    Select plot
                                </MenuItem>
                            )}
                            {availablePlots.map((plot, index) => (
                                <MenuItem key={index} value={plot}>
                                    {plot}
                                </MenuItem>
                            ))}
                        </TextField>
                    }
                    action={
                        <Box sx={{ display: 'flex', alignItems: 'center', flexShrink: 0, ml: 1 }}>
                            <BtnTooltip text={isLocked ? "Unlock dragging" : "Lock dragging"}>
                                <IconButton onClick={() => setIsLocked(!isLocked)} size="small">
                                    {isLocked ? <Lock /> : <LockOpen />}
                                </IconButton>
                            </BtnTooltip>
                            <BtnTooltip text="Add another card">
                                <IconButton
                                    color="secondary"
                                    size="small"
                                    onClick={addAnotherCard}
                                >
                                    <ContentCopy />
                                </IconButton>
                            </BtnTooltip>
                            <BtnTooltip text="Close card">
                                <IconButton onClick={closeThisCard} size="small">
                                    <Close />
                                </IconButton>
                            </BtnTooltip>
                        </Box>
                    }
                />
                <CardContent sx={{ flexGrow: 1, p: 2, overflow: 'hidden' }}>
                    {plotType === "plotly" && plotData && plotLayout ? (
                        <Plot
                            data={plotData}
                            layout={{
                                ...plotLayout,
                                autosize: true, // Keep autosize here
                                // Explicitly define margins to reduce whitespace
                                margin: {
                                    l: 20, // left margin
                                    r: 20, // right margin
                                    b: 20, // bottom margin
                                    t: 20, // top margin
                                    pad: 0 // padding between plot area and legend/color bar
                                },
                                // You might still want automargin for axes if labels can be long
                                xaxis: { automargin: true },
                                yaxis: { automargin: true },
                            }}
                            useResizeHandler={true}
                            style={{ width: "100%", height: "100%" }}
                            onClick={onPlotClick}
                            onSelected={onPlotSelected}
                            onDeselect={onPlotDeselect}
                            className="responsive-plotly-plot"
                        />
                    ) : plotType === "zndraw.Figure" && plotData ? (
                        <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%', width: '100%' }}>
                            <img
                                src={`data:image/png;base64, ${plotData}`}
                                alt="plot"
                                style={{ maxWidth: "100%", maxHeight: "100%", objectFit: "contain" }}
                            />
                        </Box>
                    ) : (
                        <Typography variant="h6" color="text.secondary" sx={{ m: 3, textAlign: 'center' }}>
                            {selectedOption ? "Loading plot data..." : "No plot selected or available."}
                        </Typography>
                    )}
                </CardContent>
            </Card>
        </Rnd>
    );
};