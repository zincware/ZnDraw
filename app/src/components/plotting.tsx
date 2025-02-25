import { useEffect, useRef, useState } from "react";
// Plottly interface
import { Button, Card, Form } from "react-bootstrap";
import { FaLock, FaLockOpen } from "react-icons/fa";
import { IoDuplicate } from "react-icons/io5";
import Plot from "react-plotly.js";
import { decodeTypedArraySpec } from "plotly.js/src/lib/array.js";
import { Rnd, type RndResizeCallback } from "react-rnd";
import * as znsocket from "znsocket";
import { client } from "../socket";
import { BtnTooltip } from "./tooltips";
import type { IndicesState } from "./utils";

interface PlottingProps {
	setStep: (step: number) => void;
	step: number;
	setSelectedFrames: (selectedFrames: IndicesState) => void;
	addPlotsWindow: number;
	setSelectedIds: (selectedIds: Set<number>) => void;
	token: string;
	updatedPlotsList: string[];
}

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
					prevCards.length > 0 ? prevCards[prevCards.length - 1] + 1 : 0;
				return [...prevCards, newCardIndex];
			});
		}
	}, [addPlotsWindow]);

	// test if a key in updatedPlotsList is not in visible plots, then create a new card
	useEffect(() => {
		updatedPlotsList.forEach((plot) => {
			if (!Object.values(visiblePlots).includes(plot)) {
				setDisplayedCards((prevCards) => {
					const newCardIndex =
						prevCards.length > 0 ? prevCards[prevCards.length - 1] + 1 : 0;
					// use visiblePlots to set the default value of the plot
					setVisiblePlots((prev: any) => {
						return { ...prev, [newCardIndex]: plot };
					});
					return [...prevCards, newCardIndex];
				});
			}
		});
	}, [updatedPlotsList]);

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
}: any) => {
	const [conInterface, setConInterface]: any =
		useState<znsocket.Dict>(undefined);
	const [availablePlots, setAvailablePlots] = useState<string[]>([]);
	const [selectedOption, setSelectedOption] = useState<string>("");
	const [rawPlotData, setRawPlotData] = useState<{ [key: string]: any }>(
		undefined,
	);
	const [plotData, setPlotData] = useState<{ [key: string]: any }>(undefined);
	const [plotType, setPlotType] = useState<string>("");
	const [plotLayout, setPlotLayout] = useState<{ [key: string]: any }>(
		undefined,
	);
	const [allowDrag, setAllowDrag] = useState<boolean>(true);
	const selectFormRef = useRef<HTMLSelectElement>(null);
	const cardRef = useRef<HTMLSelectElement>(null);

	// when created, iterate through availablePlots and set the first one as selectedOption that is not already in visiblePlots
	useEffect(() => {
		if (selectedOption !== "") {
			return;
		}
		for (let i = 0; i < availablePlots.length; i++) {
			if (!Object.values(visiblePlots).includes(availablePlots[i])) {
				setSelectedOption(availablePlots[i]);
				break;
			}
		}
	}, [availablePlots, selectedOption]);

	useEffect(() => {
		// check if identifier is in visiblePlots, if so, set selectedOption to visiblePlots[identifier]
		if (visiblePlots[identifier]) {
			setSelectedOption(visiblePlots[identifier]);
		}
	}, [identifier]);

	useEffect(() => {
		const con = new znsocket.Dict({
			client: client,
			key: `room:${token}:figures`,
		});

		con.onRefresh(async (x: any) => {
			con.keys().then((keys: any) => {
				setAvailablePlots(keys);
			});
		});

		con.keys().then((keys: any) => {
			setAvailablePlots(keys);
		});
		setConInterface(con);

		return () => {
			con.offRefresh();
		};
	}, [token]);

	useEffect(() => {
		if (conInterface === undefined) {
			return;
		}
		conInterface.get(selectedOption).then((data: any) => {
			if (data === null) {
				return;
			}
			if (data._type === "plotly.graph_objs.Figure") {
				setPlotType("plotly");
				setRawPlotData(JSON.parse(data.value).data);
				setPlotLayout(JSON.parse(data.value).layout);
			} else if (data._type === "zndraw.Figure") {
				setPlotType("zndraw.Figure");
				setPlotData(data.value.base64);
			}
		});
	}, [conInterface, selectedOption]);

	useEffect(() => {
		setVisiblePlots((prev: any) => {
			return { ...prev, [identifier]: selectedOption };
		});
	}, [selectedOption]);

	// update the actual plot if it is in updatedPlotsList
	useEffect(() => {
		if (updatedPlotsList.includes(selectedOption)) {
			conInterface.get(selectedOption).then((data: any) => {
				if (data === null) {
					return;
				}
				setRawPlotData(JSON.parse(data.value).data);
				setPlotLayout(JSON.parse(data.value).layout);
			});
		}
	}, [updatedPlotsList]);

	useEffect(() => {
		if (rawPlotData && plotType === "plotly") {
			const markerList: [number, number, string][] = [];

			// Function to convert fields containing bdata and dtype
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

			// Process each data item
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
						if (customdata[0] === step) {
							let xPosition = convertedItem.x[index];
							let yPosition = convertedItem.y[index];

							let color = "red";
							if (convertedItem.line?.color) {
								color = convertedItem.line.color;
							}

							markerList.push([xPosition, yPosition, color]);
						}
					});
				}

				return convertedItem;
			});

			// Add marker data
			updatedPlotData.push({
				type: "scatter",
				mode: "markers",
				name: "Step",
				showlegend: false,
				x: markerList.map((marker) => marker[0]),
				y: markerList.map((marker) => marker[1]),
				marker: {
					color: markerList.map((marker) => marker[2]),
					size: 10,
					symbol: "circle",
					line: {
						color: "black",
						width: 2,
					},
				},
			});

			setPlotData(updatedPlotData);
		}
	}, [rawPlotData, step, plotType]);

	const handleSelectChange = (event: React.ChangeEvent<HTMLSelectElement>) => {
		setSelectedOption(event.target.value);
	};

	const closeThisCard = () => {
		setDisplayedCards((prevCards) =>
			prevCards.filter((card) => card !== identifier),
		);
		// also remove from visiblePlots
		setVisiblePlots((prev: any) => {
			const copy = { ...prev };
			delete copy[identifier];
			return copy;
		});
	};

	const onResize: RndResizeCallback = () => {
		if (cardRef.current) {
			setPlotLayout((prev) => {
				if (prev) {
					return {
						...prev,
						width: cardRef.current.clientWidth - 20,
						height: cardRef.current.clientHeight - 60,
					};
				}
				return prev;
			});
		}
	};

	const onPlotClick = ({ points }: { points: any[] }) => {
		if (points[0]?.customdata[0]) {
			setStep(points[0].customdata[0]);
		}
		if (points[0]?.customdata[1]) {
			setSelectedIds(new Set([points[0].customdata[1]]));
		}
	};

	const addAnotherCard = () => {
		setDisplayedCards((prevCards) => {
			const newCardIndex =
				prevCards.length > 0 ? prevCards[prevCards.length - 1] + 1 : 0;
			return [...prevCards, newCardIndex];
		});
	};

	const onPlotSelected = (event: any) => {
		if (!event || !event.points) {
			return;
		}
		if (event.points.length === 0) {
			// This is triggered once the plot is re-rendered. We want to keep the selection here.
			return;
		}
		const selectedFrames = event.points.map((point: any) =>
			point.customdata ? point.customdata[0] : point.pointIndex,
		);
		// for all points.customdata[0] == step collect the points.customdata[1] and set selectedIds if customdata[1] is available
		const selectedIds = new Set<number>(
			event.points
				.filter((point: any) => point.customdata?.[1])
				.map((point: any) => point.customdata[1]),
		);
		if (selectedIds.size > 0) {
			setSelectedIds(selectedIds);
		}

		setSelectedFrames({
			active: true,
			indices: new Set(selectedFrames),
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
			onResize={onResize}
			disableDragging={!allowDrag}
			dragGrid={[50, 50]}
			resizeGrid={[50, 50]}
		>
			<Card
				style={{
					padding: 0,
					width: "100%",
					height: "100%",
				}}
				ref={cardRef}
			>
				<Card.Header
					className="d-flex justify-content-between align-items-center flex-nowrap"
					style={{ height: 50 }}
					onPointerEnter={() => setAllowDrag(true)}
					onPointerLeave={() => setAllowDrag(false)}
				>
					<Form.Select
						onChange={handleSelectChange}
						value={selectedOption} // https://github.com/react-bootstrap/react-bootstrap/issues/2091
						ref={selectFormRef}
					>
						{selectedOption === "" && (
							<option value="" disabled>
								Select plot
							</option>
						)}
						{availablePlots.map((plot, index) => (
							<option key={index} value={plot}>
								{plot}
							</option>
						))}
					</Form.Select>
					<BtnTooltip text="Add another card">
						<Button
							variant="tertiary"
							className="mx-2 btn btn-outline-secondary"
							onClick={addAnotherCard}
						>
							<IoDuplicate />
						</Button>
					</BtnTooltip>
					<Button variant="close" className="mx-2" onClick={closeThisCard} />
				</Card.Header>
				<Card.Body style={{ padding: 0 }}>
					{plotType === "plotly" && (
						<Plot
							data={plotData}
							layout={plotLayout}
							onClick={onPlotClick}
							onSelected={onPlotSelected}
							onDeselect={onPlotDeselect}
						/>
					)}
					{plotType === "zndraw.Figure" && (
						<img
							src={`data:image/png;base64, ${plotData}`}
							alt="plot"
							className="img-fluid"
						/>
					)}
					{plotType === "" && (
						<h3 className="text-secondary m-3">No data available</h3>
					)}
				</Card.Body>
			</Card>
		</Rnd>
	);
};
