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
  let [conInterface, setConInterface]: any = useState<znsocket.Dict>(undefined);
  let [availablePlots, setAvailablePlots] = useState<string[]>([]);
  let [selectedOption, setSelectedOption] = useState<string>("");
  let [rawPlotData, setRawPlotData] = useState<{ [key: string]: any }>(
    undefined,
  );
  let [plotData, setPlotData] = useState<{ [key: string]: any }>(undefined);
  let [plotType, setPlotType] = useState<string>("");
  let [plotLayout, setPlotLayout] = useState<{ [key: string]: any }>(undefined);
  let [plotHover, setPlotHover] = useState<boolean>(false);
  const [allowDrag, setAllowDrag] = useState<boolean>(true);
  let selectFormRef = useRef<HTMLSelectElement>(null);
  let cardRef = useRef<HTMLSelectElement>(null);

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
    if (plotHover) {
      setAllowDrag(false);
    } else {
      const debounceTimeout = setTimeout(() => {
        setAllowDrag(true);
      }, 1000);
      return () => clearTimeout(debounceTimeout);
    }
  }, [plotHover]);

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
    conInterface.get(selectedOption).then((data: any) => {
      if (data === null) {
        return;
      }
      if (data["_type"] === "plotly.graph_objs.Figure") {
        setPlotType("plotly");
        setRawPlotData(JSON.parse(data["value"]).data);
        setPlotLayout(JSON.parse(data["value"]).layout);
      } else if (data["_type"] === "zndraw.Figure") {
        setPlotType("zndraw.Figure");
        setPlotData(data["value"]["base64"]);
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
        setRawPlotData(JSON.parse(data["value"]).data);
        setPlotLayout(JSON.parse(data["value"]).layout);
      });
    }
  }, [updatedPlotsList]);

  useEffect(() => {
    if (rawPlotData && plotType === "plotly") {
      const markerList: [number, number, string][] = [];

      // Add markers at the matching step in the data
      rawPlotData.forEach((dataItem) => {
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

      const plotDataCopy = JSON.parse(JSON.stringify(rawPlotData));

      // Add the markers to the data array
      plotDataCopy.push({
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
      setPlotData(plotDataCopy);
    }
  }, [rawPlotData, step, plotType]); // does this self-trigger? If so use raw

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
          {plotType == "plotly" && (
            <Plot
              data={plotData}
              layout={plotLayout}
              onHover={() => setPlotHover(true)}
              onSelecting={() => setPlotHover(true)}
              onBeforeHover={() => setPlotHover(true)}
              onUnhover={() => setPlotHover(false)}
              onClick={onPlotClick}
              onSelected={onPlotSelected}
              onDeselect={onPlotDeselect}
            />
          )}
          {plotType == "zndraw.Figure" && (
            <img
              src={`data:image/png;base64, ${plotData}`}
              alt="plot"
              className="img-fluid"
            />
          )}
          {plotType == "" && (
            <h3 className="text-secondary m-3">No data available</h3>
          )}
        </Card.Body>
      </Card>
    </Rnd>
  );
};
