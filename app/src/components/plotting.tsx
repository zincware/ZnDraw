// Plottly interface
import { Button, Card, Form } from "react-bootstrap";
import { useState, useRef, useEffect } from "react";
import { socket } from "../socket";
import { Rnd, RndResizeCallback } from "react-rnd";
import Plot from "react-plotly.js";
import { IoDuplicate } from "react-icons/io5";

interface PlottingProps {
  setStep: (step: number) => void;
  setSelectedFrames: (selectedFrames: Set<number>) => void;
  addPlotsWindow: number;
}

export const Plotting = ({ setStep, setSelectedFrames, addPlotsWindow }: PlottingProps) => {
  const [availablePlots, setAvailablePlots] = useState<string[]>([]);
  const [plotData, setPlotData] = useState<{ [key: string]: any }>({});
  const [displayedCards, setDisplayedCards] = useState<number[]>([]);

  useEffect(() => {
    if (addPlotsWindow > 0) {
      setDisplayedCards((prevCards) => {
        const newCardIndex =
          prevCards.length > 0 ? prevCards[prevCards.length - 1] + 1 : 0;
        return [...prevCards, newCardIndex];
      });
    }
  }, [addPlotsWindow]);

  // on analysis:figure:refresh add another card
  useEffect(() => {
    const handleFigureRefresh = () => {
      setDisplayedCards((prevCards) => {
        const newCardIndex =
          prevCards.length > 0 ? prevCards[prevCards.length - 1] + 1 : 0;
        return [...prevCards, newCardIndex];
      });
    };
    socket.on("analysis:figure:refresh", handleFigureRefresh);
    return () => {
      socket.off("analysis:figure:refresh", handleFigureRefresh);
    };
  }, []); // Removed displayedCards from dependencies

  useEffect(() => {
    availablePlots.forEach((plot) => {
      socket.emit("analysis:figure:get", plot, (data: any) => {
        setPlotData((prevData) => ({
          ...prevData,
          [plot]: JSON.parse(data),
        }));
      });
    });
  }, [availablePlots]);

  return (
    <>
      {displayedCards.map((cardIndex) => (
        <PlotsCard
          key={cardIndex}
          identifier={cardIndex}
          availablePlots={availablePlots}
          setAvailablePlots={setAvailablePlots}
          plotData={plotData}
          setDisplayedCards={setDisplayedCards}
          setStep={setStep}
          setSelectedFrames={setSelectedFrames}
        />
      ))}
    </>
  );
};

interface PlotsCardProps {
  identifier: number;
  availablePlots: string[];
  setAvailablePlots: (availablePlots: string[]) => void;
  plotData: { [key: string]: any };
  setDisplayedCards: (displayedCards: number[]) => void;
  setStep: (step: number) => void;
  setSelectedFrames: (selectedFrames: Set<number>) => void;
}

const PlotsCard = ({
  identifier,
  availablePlots,
  setAvailablePlots,
  plotData,
  setDisplayedCards,
  setStep,
  setSelectedFrames,
}: PlotsCardProps) => {
  const cardRef = useRef<any>(null);
  const [selectedOption, setSelectedOption] = useState<string>("");
  const [allowDrag, setAllowDrag] = useState<boolean>(true);
  const [plotLayout, setPlotLayout] = useState<any>({});

  const handleSelectChange = (event: React.ChangeEvent<HTMLSelectElement>) => {
    setSelectedOption(event.target.value);
  };

  // set initial data if availablePlots is not empty
  useEffect(() => {
    socket.emit("analysis:figure:keys", (data: string[]) => {
      setAvailablePlots(data);
    });
  }, []);

  // once plot data updates and selectedOption == "" set selectedOption to first available plot
  // TODO: this part is still very buggy!
  useEffect(() => {
    if (availablePlots.length > 0 && selectedOption === "") {
      setSelectedOption(availablePlots[0]);
    }
  }, [availablePlots, selectedOption]);

  useEffect(() => {
    if (plotData[selectedOption]) {
      setPlotLayout(plotData[selectedOption].layout);
    }
  }, [plotData, selectedOption]);

  const onPlotClick = ({ points }: { points: any[] }) => {
    if (points[0]?.customdata) {
      setStep(points[0].customdata[0]);
    } else {
      setStep(points[0].pointIndex);
    }
  };

  const onPlotSelected = (event: any) => {
    const selectedFrames = event.points.map((point: any) =>
      point.customdata ? point.customdata[0] : point.pointIndex,
    );
    setSelectedFrames(new Set(selectedFrames));
  };

  const handleSelectClick = () => {
    socket.emit("analysis:figure:keys", (data: string[]) => {
      setAvailablePlots(data);
    });
  };

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
          <Button
            variant="tertiary"
            className="mx-2 btn btn-outline-secondary"
            onClick={addAnotherCard}
          >
            <IoDuplicate />
          </Button>
          <Button variant="close" className="mx-2" onClick={closeThisCard} />
        </Card.Header>
        <Card.Body style={{ padding: 0 }}>
          {plotData[selectedOption] ? (
            <Plot
              data={plotData[selectedOption].data}
              frames={plotData[selectedOption].frames}
              config={plotData[selectedOption].config}
              layout={plotLayout}
              onHover={() => setAllowDrag(false)}
              onUnhover={() => setAllowDrag(true)}
              onClick={onPlotClick}
              onSelected={onPlotSelected}
              onDeselect={() => setSelectedFrames(new Set())}
            />
          ) : (
            <h3 className="text-secondary m-3">No data available</h3>
          )}
        </Card.Body>
      </Card>
    </Rnd>
  );
};
