// Plottly interface
import { Button, Card, Form } from "react-bootstrap";
import { useState, useRef, useEffect } from "react";
import { socket } from "../socket";
import { Rnd, RndResizeCallback } from "react-rnd";
import Plot from "react-plotly.js";
import { IoDuplicate } from "react-icons/io5";

export const Plotting = () => {

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
        />
      ))}

      {/* <PlotsCard 
      availablePlots={availablePlots}
      setAvailablePlots={setAvailablePlots}
      plotData={plotData}
      /> */}
      {/* <PlotsCard /> */}
      {/* <PlotsCard />
      <PlotsCard /> */}
    </>
  );
};

interface handleFigureDataProps {
  setPlotData: (data: any) => void;
  setShowPlotsCard: (value: boolean) => void;
}

const handleFigureData = ({
  setPlotData,
  setShowPlotsCard,
}: handleFigureDataProps) => {
  useEffect(() => {
    const handleFigure = (data) => {
      try {
        const parsedData = JSON.parse(data);
        setPlotData(parsedData);
        setShowPlotsCard(true);
      } catch (error) {
        console.error("Error parsing JSON data: ", error);
      }
    };

    socket.on("analysis:figure:set", handleFigure);
    return () => {
      socket.off("analysis:figure:set", handleFigure);
    };
  }, []);
};

const PlotsCard = ({identifier, availablePlots, setAvailablePlots, plotData, setDisplayedCards}: any) => {
  const cardRef = useRef<any>(null);
  const [selectedOption, setSelectedOption] = useState<string>("1");
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

  const handleSelectClick = () => {
    // call 'analysis:figure:keys' to get the available plots
    socket.emit("analysis:figure:keys", (data: string[]) => {
      setAvailablePlots(data);
    });
  }

  const closeThisCard = () => {
    // remove this card from the displayed cards
    setDisplayedCards((prevCards: number[]) => {
      return prevCards.filter((card: number) => card !== identifier);
    });
  }

  const addAnotherCard = () => {
    // add another card to the displayed cards
    setDisplayedCards((prevCards: number[]) => {
      return [...prevCards, prevCards[prevCards.length - 1] + 1];
    });
  }

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
          {/* <Card.Title className="p-2">Analysis</Card.Title> */}
          <Form.Select onChange={handleSelectChange} onClick={handleSelectClick}>
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
          <Button variant="close" className="mx-2" onClick={closeThisCard}/>
        </Card.Header>
        <Card.Body style={{ padding: 0 }}>
        {plotData[selectedOption] && (
              <Plot
                data={plotData[selectedOption].data}
                frames={plotData[selectedOption].frames}
                config={plotData[selectedOption].config}
                layout={plotLayout} // todo: merge
                onBeforeHover={() => setAllowDrag(false)}
                onUnhover={() => setAllowDrag(true)}
              />
            )}
          {/* {selectedOption === "1" && (
            <Plot
              data={[
                {
                  x: [1, 2, 3],
                  y: [2.5, 6.5, 3.5],
                  type: "scatter",
                  mode: "lines+markers",
                  marker: { color: "red" },
                },
                { type: "bar", x: [1, 2, 3], y: [2, 5, 3] },
              ]}
              layout={plotLayout}
              onBeforeHover={() => setAllowDrag(false)}
              onUnhover={() => setAllowDrag(true)}
            />
          )} */}
        </Card.Body>
      </Card>
    </Rnd>
  );
};

//   const PlotsCard = ({
//     plotData,
//     setPlotData,
//     colorMode,
//     showPlotsCard,
//     setShowPlotsCard,
//     setStep,
//   }: {
//     plotData: any;
//     setPlotData: any;
//     colorMode: string;
//     showPlotsCard: boolean;
//     setShowPlotsCard: any;
//     setStep: any;
//   }) => {
//     const [plotStyle, setPlotStyle] = useState<any>({
//       width: "100%",
//       height: "100%",
//     });
//     const [renderKey, setRenderKey] = useState<number>(0);
//     const cardRef = useRef<any>(null);

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

//     const onResize = () => {
//       if (cardRef.current) {
//         const newPlotData = { ...plotData };
//         newPlotData.layout.width = cardRef.current.clientWidth;
//         newPlotData.layout.height = cardRef.current.clientHeight - 50;
//         setPlotData(newPlotData);
//         setRenderKey((prevKey) => prevKey + 1);
//       }
//     };

//     const onPlotClick = ({
//       event,
//       points,
//     }: {
//       event: MouseEvent;
//       points: any[];
//     }) => {
//       if (points[0].customdata) {
//         setStep(points[0].customdata[0]);
//       } else {
//         setStep(points[0].pointIndex);
//       }
//     };

//     return (
//       <Rnd
//         default={{
//           x: 100,
//           y: -100,
//           width: 400,
//           height: 400,
//         }}
//         minHeight={200}
//         minWidth={220}
//         style={{
//           zIndex: 1000,
//           padding: 0,
//           margin: 0,
//           display: showPlotsCard ? "block" : "none",
//         }}
//         onResize={onResize}
//       >
//         <Card
//           style={{
//             margin: 0,
//             padding: 0,
//             width: "100%",
//             height: "100%",
//           }}
//           ref={cardRef}
//         >
//           <Card.Header className="d-flex justify-content-between align-items-center">
//             <Card.Title>Analysis Figure</Card.Title>
//             <Button variant="close" onClick={() => setShowPlotsCard(false)} />
//           </Card.Header>
//           <Card.Body>
//             {plotData.data.length > 0 && (
//               <Plot
//                 key={renderKey}
//                 data={plotData.data}
//                 layout={plotData.layout}
//                 frames={plotData.frames}
//                 config={plotData.config}
//                 style={plotStyle}
//                 onClick={onPlotClick}
//               />
//             )}
//           </Card.Body>
//         </Card>
//       </Rnd>
//     );
//   };
