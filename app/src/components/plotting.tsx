// Plottly interface
import { Button, Card, Form } from "react-bootstrap";
import { useState, useRef, useEffect } from "react";
import { client } from "../socket";
import * as znsocket from "znsocket";
import { Rnd, RndResizeCallback } from "react-rnd";
import Plot from "react-plotly.js";
import { IoDuplicate } from "react-icons/io5";
import { FaLock, FaLockOpen } from "react-icons/fa";
import { BtnTooltip } from "./tooltips";
import { IndicesState } from "./utils";

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
  // const [availablePlots, setAvailablePlots] = useState<string[]>([]);
  // const [plotData, setPlotData] = useState<{ [key: string]: any }>({});
  const [displayedCards, setDisplayedCards] = useState<number[]>([]);
  const [visiblePlots, setVisiblePlots] = useState<{ [key: number]: string }>({});

  // TODO: if the displayedCards is closed, visiblePlots is not beeing removed properly

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
      // TODO: we need to set a "default plot" to be displayed when a new plot is added
      if (!Object.values(visiblePlots).includes(plot)) {
        setDisplayedCards((prevCards) => {
          const newCardIndex =
            prevCards.length > 0 ? prevCards[prevCards.length - 1] + 1 : 0;
          // use visiblePlots to set the default value of the plot
          setVisiblePlots((prev: any) => {
            return {...prev, [newCardIndex]: plot};
          });
          return [...prevCards, newCardIndex];
        });
      }
    });
  }, [updatedPlotsList]);

  // setAddPlotsWindow((prev: number) => prev + 1)}

  // on analysis:figure:refresh add another card
  // useEffect(() => {
  //   const handleFigureRefresh = () => {
  //     setDisplayedCards((prevCards) => {
  //       const newCardIndex =
  //         prevCards.length > 0 ? prevCards[prevCards.length - 1] + 1 : 0;
  //       return [...prevCards, newCardIndex];
  //     });
  //   };
  //   socket.on("analysis:figure:refresh", handleFigureRefresh);
  //   return () => {
  //     socket.off("analysis:figure:refresh", handleFigureRefresh);
  //   };
  // }, []); // Removed displayedCards from dependencies

  // useEffect(() => {
  //   availablePlots.forEach((plot) => {
  //     socket.emit("analysis:figure:get", plot, (data: any) => {
  //       setPlotData((prevData) => ({
  //         ...prevData,
  //         [plot]: JSON.parse(data),
  //       }));
  //     });
  //   });
  // }, [availablePlots]);

  return (
    <>
    {displayedCards.map((cardIndex) => (
      <PlotsCard2 key={cardIndex} identifier={cardIndex} updatedPlotsList={updatedPlotsList} token={token} setVisiblePlots={setVisiblePlots} setDisplayedCards={setDisplayedCards} visiblePlots={visiblePlots} setStep={setStep} setSelectedIds={setSelectedIds} setSelectedFrames={setSelectedFrames} step={step}/>
    ))}

      {/* {displayedCards.map((cardIndex) => (
        <PlotsCard
          key={cardIndex}
          identifier={cardIndex}
          availablePlots={availablePlots}
          setAvailablePlots={setAvailablePlots}
          plotData={plotData}
          // current plot data
          setDisplayedCards={setDisplayedCards}
          setStep={setStep}
          setSelectedFrames={setSelectedFrames}
          setSelectedIds={setSelectedIds}
          step={step}
        />
      ))} */}
    </>
  );
};


const PlotsCard2 = ({updatedPlotsList, token, setVisiblePlots, identifier, setDisplayedCards, visiblePlots, setStep, setSelectedIds, setSelectedFrames, step}: any) => {
  let [conInterface, setConInterface]: any = useState(undefined);
  let [availablePlots, setAvailablePlots] = useState<string[]>([]);
  let [selectedOption, setSelectedOption] = useState<string>("");
  let [plotData, setPlotData] = useState<{ [key: string]: any }>(undefined);
  const [allowDrag, setAllowDrag] = useState<boolean>(true);
  let selectFormRef = useRef<HTMLSelectElement>(null);

  useEffect(() => {
    // check if identifier is in visiblePlots, if so, set selectedOption to visiblePlots[identifier]
    if (visiblePlots[identifier]) {
      setSelectedOption(visiblePlots[identifier]);
    }
  }, [identifier]);

  useEffect(() => {
    const con = new znsocket.Dict({
      client: client,
      key: "room:" + token + ":figures",
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
    conInterface.getitem(selectedOption).then((data: any) => {
      if (data === null) {
        return;
      }
      setPlotData(JSON.parse(data["value"]));
    });
  }, [conInterface, selectedOption]);

  useEffect(() => {
    setVisiblePlots((prev: any) => {
      return {...prev, [identifier]: selectedOption};
    });
  }, [selectedOption]);
  
  // update the actuall plot if it is in updatedPlotsList
  useEffect(() => {
    if (updatedPlotsList.includes(selectedOption)) {
      conInterface.getitem(selectedOption).then((data: any) => {
        if (data === null) {
          return;
        }
        setPlotData(JSON.parse(data["value"]));
      });
    }
  }, [updatedPlotsList]);

  // useEffect(() => {
  //   if (plotData[selectedOption]) {
  //     if (plotData[selectedOption].layout) {
  //       // Deep copy the layout and data to prevent mutating the original
  //       const layout = JSON.parse(
  //         JSON.stringify(plotData[selectedOption].layout),
  //       );
  //       const data = JSON.parse(JSON.stringify(plotData[selectedOption].data));

  //       const markerList: [number, number, string][] = [];

  //       // Add markers at the matching step in the data
  //       plotData[selectedOption].data.forEach((dataItem) => {
  //         if (dataItem.customdata) {
  //           dataItem.customdata.forEach((customdata, index) => {
  //             // Check if customdata[0] matches the step
  //             if (customdata[0] === step) {
  //               const xPosition = dataItem.x[index];
  //               const yPosition = dataItem.y[index];
  //               // check if dataItem.line.color is available
  //               let color = "red";
  //               if (dataItem.line) {
  //                 if (dataItem.line.color) {
  //                   color = dataItem.line.color;
  //                 }
  //               }
  //               markerList.push([xPosition, yPosition, color]);
  //             }
  //           });
  //         }
  //       });

  //       // Add the markers to the data array
  //       data.push({
  //         type: "scatter",
  //         mode: "markers",
  //         name: "Step",
  //         showlegend: false,
  //         x: markerList.map((marker) => marker[0]),
  //         y: markerList.map((marker) => marker[1]),
  //         marker: {
  //           color: markerList.map((marker) => marker[2]),
  //           size: 10,
  //           symbol: "circle",
  //           line: {
  //             color: "black",
  //             width: 2,
  //           },
  //         },
  //       });

  //       // Set the updated layout and data
  //       setPlotLayout(layout);
  //       setActualPlotData(data);
  //     }
  //   } else {
  //     setActualPlotData(null);
  //   }
  // }, [plotData, selectedOption, step]);

  const handleSelectChange = (event: React.ChangeEvent<HTMLSelectElement>) => {
    setSelectedOption(event.target.value);
  };

  const closeThisCard = () => {
    setDisplayedCards((prevCards) =>
      prevCards.filter((card) => card !== identifier),
    );
    // also remove from visiblePlots
    setVisiblePlots((prev: any) => {
      const copy = {...prev};
      delete copy[identifier];
      return copy;
    });
  };

  const onPlotClick = ({ points }: { points: any[] }) => {
    if (points[0]?.customdata) {
      setStep(points[0].customdata[0]);
      if (points[0].customdata[1]) {
        setSelectedIds(new Set([points[0].customdata[1]]));
      }
    }
  };

  const onPlotSelected = (event: any) => {
    setAllowDrag(false);
    const selectedFrames = event.points.map((point: any) =>
      point.customdata ? point.customdata[0] : point.pointIndex,
    );
    // for all points.customdata[0] == step collect the points.customdata[1] and set selectedIds if customdata[1] is available
    const selectedIds = new Set<number>(
      event.points
        .filter((point: any) => point.customdata && point.customdata[1])
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
      // onResize={onResize}
      // disableDragging={!allowDrag}
      dragGrid={[50, 50]}
      resizeGrid={[50, 50]}
    >
      <Card
        style={{
          padding: 0,
          width: "100%",
          height: "100%",
        }}
        // ref={cardRef}
      >
        <Card.Header
          className="d-flex justify-content-between align-items-center flex-nowrap"
          style={{ height: 50 }}
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
          {/* <BtnTooltip
            text={allowDrag ? "Lock card movement" : "Unlock card movement"}
          >
            <Button
              variant="outline-secondary"
              className="mx-1"
              onClick={() => setAllowDrag(!allowDrag)}
            >
              {allowDrag ? <FaLockOpen /> : <FaLock />}
            </Button>
          </BtnTooltip> */}
          <BtnTooltip text="Add another card">
            <Button
              variant="tertiary"
              className="mx-2 btn btn-outline-secondary"
              // onClick={addAnotherCard}
            >
              <IoDuplicate />
            </Button>
          </BtnTooltip>
          <Button variant="close" className="mx-2" onClick={closeThisCard} />
        </Card.Header>
        <Card.Body style={{ padding: 0 }}>
          {plotData ? (
            <Plot
              data={plotData.data}
              // frames={plotData[selectedOption].frames}
              // config={plotData[selectedOption].config}
              // layout={plotLayout}
              onHover={() => setAllowDrag(false)}
              onUnhover={() => setAllowDrag(true)}
              onClick={onPlotClick}
              onSelected={onPlotSelected}
              onDeselect={onPlotDeselect}
            />
          ) : (
            <h3 className="text-secondary m-3">No data available</h3>
          )}
        </Card.Body>
      </Card>
    </Rnd>
  );
};


interface PlotsCardProps {
  identifier: number;
  availablePlots: string[];
  // setAvailablePlots: (availablePlots: string[]) => void;
  plotData: { [key: string]: any };
  setDisplayedCards: (displayedCards: number[]) => void;
  setStep: (step: number) => void;
  step: number;
  setSelectedFrames: (selectedFrames: IndicesState) => void;
  setSelectedIds: (selectedIds: Set<number>) => void;
}

const PlotsCard = ({
  identifier,
  availablePlots,
  // setAvailablePlots,
  plotData,
  setDisplayedCards,
  setStep,
  step,
  setSelectedFrames,
  setSelectedIds,
}: PlotsCardProps) => {
  const cardRef = useRef<any>(null);
  const [selectedOption, setSelectedOption] = useState<string>("");
  const [allowDrag, setAllowDrag] = useState<boolean>(true);
  const [plotLayout, setPlotLayout] = useState<any>({});
  const [actualPlotData, setActualPlotData] = useState<any>(null);

  const handleSelectChange = (event: React.ChangeEvent<HTMLSelectElement>) => {
    setSelectedOption(event.target.value);
  };


  useEffect(() => {
    if (plotData[selectedOption]) {
      if (plotData[selectedOption].layout) {
        // Deep copy the layout and data to prevent mutating the original
        const layout = JSON.parse(
          JSON.stringify(plotData[selectedOption].layout),
        );
        const data = JSON.parse(JSON.stringify(plotData[selectedOption].data));

        const markerList: [number, number, string][] = [];

        // Add markers at the matching step in the data
        plotData[selectedOption].data.forEach((dataItem) => {
          if (dataItem.customdata) {
            dataItem.customdata.forEach((customdata, index) => {
              // Check if customdata[0] matches the step
              if (customdata[0] === step) {
                const xPosition = dataItem.x[index];
                const yPosition = dataItem.y[index];
                // check if dataItem.line.color is available
                let color = "red";
                if (dataItem.line) {
                  if (dataItem.line.color) {
                    color = dataItem.line.color;
                  }
                }
                markerList.push([xPosition, yPosition, color]);
              }
            });
          }
        });

        // Add the markers to the data array
        data.push({
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

        // Set the updated layout and data
        setPlotLayout(layout);
        setActualPlotData(data);
      }
    } else {
      setActualPlotData(null);
    }
  }, [plotData, selectedOption, step]);

  const onPlotClick = ({ points }: { points: any[] }) => {
    if (points[0]?.customdata) {
      setStep(points[0].customdata[0]);
      if (points[0].customdata[1]) {
        setSelectedIds(new Set([points[0].customdata[1]]));
      }
    }
  };

  const onPlotSelected = (event: any) => {
    setAllowDrag(false);
    const selectedFrames = event.points.map((point: any) =>
      point.customdata ? point.customdata[0] : point.pointIndex,
    );
    // for all points.customdata[0] == step collect the points.customdata[1] and set selectedIds if customdata[1] is available
    const selectedIds = new Set<number>(
      event.points
        .filter((point: any) => point.customdata && point.customdata[1])
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

  // const handleSelectClick = () => {
  //   socket.emit("analysis:figure:keys", (data: string[]) => {
  //     setAvailablePlots(data);
  //   });
  // };

  const closeThisCard = () => {
    setDisplayedCards((prevCards) =>
      prevCards.filter((card) => card !== identifier),
    );
  };

  const addAnotherCard = () => {
    setDisplayedCards((prevCards) => {
      const newCardIndex =
        prevCards.length > 0 ? prevCards[prevCards.length - 1] + 1 : 0;
      return [...prevCards, newCardIndex];
    });
  };

  const onResize: RndResizeCallback = () => {
    if (cardRef.current) {
      setPlotLayout((prev) => ({
        ...prev,
        width: cardRef.current.clientWidth - 20,
        height: cardRef.current.clientHeight - 60,
      }));
    }
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
        >
          <Form.Select
            onChange={handleSelectChange}
            onClick={handleSelectClick}
            defaultValue=""
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
          <BtnTooltip
            text={allowDrag ? "Lock card movement" : "Unlock card movement"}
          >
            <Button
              variant="outline-secondary"
              className="mx-1"
              onClick={() => setAllowDrag(!allowDrag)}
            >
              {allowDrag ? <FaLockOpen /> : <FaLock />}
            </Button>
          </BtnTooltip>
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
          {plotData[selectedOption] && actualPlotData ? (
            <Plot
              data={actualPlotData}
              frames={plotData[selectedOption].frames}
              config={plotData[selectedOption].config}
              layout={plotLayout}
              onHover={() => setAllowDrag(false)}
              // onUnhover={() => setAllowDrag(true)}
              onClick={onPlotClick}
              onSelected={onPlotSelected}
              onDeselect={onPlotDeselect}
            />
          ) : (
            <h3 className="text-secondary m-3">No data available</h3>
          )}
        </Card.Body>
      </Card>
    </Rnd>
  );
};
