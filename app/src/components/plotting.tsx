// Plottly interface
import { Button, Card, Form } from "react-bootstrap";
import { useState, useRef, useEffect } from "react";
import { socket } from "../socket";
import { Rnd, RndResizeCallback } from "react-rnd";
import Plot from "react-plotly.js";
import { IoDuplicate } from "react-icons/io5";


interface PlottingProps {
  setStep: (step: number) => void;
}

export const Plotting = ({ setStep }: PlottingProps) => {
  const [availablePlots, setAvailablePlots] = useState<string[]>([]);
  const [plotData, setPlotData] = useState<object>({}); // dict[string| dict]
  const [displayedCards, setDisplayedCards] = useState<number[]>([1]);

  useEffect(() => {
    // when availablePlots is updated, update all the plotData
    availablePlots.forEach((plot) => {
      socket.emit("analysis:figure:get", plot, (data: any) => {
        setPlotData((prevData: object) => {
          return {
            ...prevData,
            [plot]: JSON.parse(data),
          };
        });
      });
    });
    console.log("availablePlots: ", availablePlots);
    // log the keys of the plotData
    console.log("plotData keys: ", Object.keys(plotData));
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
        />
      ))}
    </>
  );
};

//     useEffect(() => {
//       if (plotData) {
//         const newPlotData = { ...plotData };
//         newPlotData.layout.paper_bgcolor =
//           colorMode === "dark" ? "rgba(0,0,0, 0)" : "rgba(255,255,255, 0)";
//         setPlotData(newPlotData);
//       }
//       const newPlotStyle = { ...plotStyle };
//       newPlotStyle.filter =
//         colorMode === "dark" ? "invert(75%) hue-rotate(180deg)" : "";
//       setPlotStyle(newPlotStyle);
//     }, [colorMode]);

interface PlotsCardProps {
  identifier: number;
  availablePlots: string[];
  setAvailablePlots: (availablePlots: string[]) => void;
  plotData: object;
  setDisplayedCards: (displayedCards: number[]) => void;
  setStep: (step: number) => void;
}

const PlotsCard = ({
  identifier,
  availablePlots,
  setAvailablePlots,
  plotData,
  setDisplayedCards,
  setStep,
}: PlotsCardProps) => {
  const cardRef = useRef<any>(null);
  const [selectedOption, setSelectedOption] = useState<string>("");
  const [allowDrag, setAllowDrag] = useState<boolean>(true);
  const [plotLayout, setPlotLayout] = useState<any>({
    width: 220 - 20,
    height: 200 - 60,
  });

  const handleSelectChange = (event: any) => {
    // update the selected plot via 'analysis:figure:get' with event.target.value
    console.log("selected option: ", event.target.value);
    setSelectedOption(event.target.value);
  };

  const onPlotClick = ({
    event,
    points,
  }: {
    event: MouseEvent;
    points: any[];
  }) => {
    if (points[0].customdata) {
      setStep(points[0].customdata[0]);
    } else {
      setStep(points[0].pointIndex);
    }
  };

  const handleSelectClick = () => {
    // call 'analysis:figure:keys' to get the available plots
    socket.emit("analysis:figure:keys", (data: string[]) => {
      setAvailablePlots(data);
    });
  };

  const closeThisCard = () => {
    // remove this card from the displayed cards
    setDisplayedCards((prevCards: number[]) => {
      return prevCards.filter((card: number) => card !== identifier);
    });
  };

  const addAnotherCard = () => {
    // add another card to the displayed cards
    setDisplayedCards((prevCards: number[]) => {
      return [...prevCards, prevCards[prevCards.length - 1] + 1];
    });
  };

  useEffect(() => {
    console.log("selected option: ", selectedOption);
    console.log("plotData: ", plotData[selectedOption]);
  }, [selectedOption]);

  const onResize: RndResizeCallback = () => {
    if (cardRef.current) {
      setPlotLayout({
        width: cardRef.current.clientWidth - 20,
        height: cardRef.current.clientHeight - 60,
        // margin of 10px on bottom and sides
      });
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
            className=" mx-2 btn btn-outline-secondary"
            onClick={addAnotherCard}
          >
            <IoDuplicate />
          </Button>
          {/* TODO: tooltip */}
          <Button variant="close" className="mx-2" onClick={closeThisCard} />
        </Card.Header>
        <Card.Body style={{ padding: 0 }}>
        {plotData[selectedOption] ? (
  <Plot
    data={plotData[selectedOption].data}
    frames={plotData[selectedOption].frames}
    config={plotData[selectedOption].config}
    layout={plotLayout} // todo: merge
    onBeforeHover={() => setAllowDrag(false)}
    onUnhover={() => setAllowDrag(true)}
    onClick={onPlotClick}
  />
) : (
  <h3 className="text-secondary m-3">No data available</h3>
)}
        </Card.Body>
      </Card>
    </Rnd>
  );
};