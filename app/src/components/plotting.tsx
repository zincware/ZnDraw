// Plottly interface
import { Button, Card, Form } from "react-bootstrap";
import { useState, useRef, useEffect } from "react";
import { socket } from "../socket";
import { Rnd, RndResizeCallback } from "react-rnd";
import Plot from "react-plotly.js";
import { all } from "three/examples/jsm/nodes/Nodes.js";
import { IoDuplicate } from "react-icons/io5";

export const Plotting = () => {
  useEffect(() => {
    console.log("Setting up plots");
  }, []);

  return (
    <>
      <PlotsCard />
      <PlotsCard />
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

const PlotsCard = () => {
  const cardRef = useRef<any>(null);
  const [selectedOption, setSelectedOption] = useState<string>("1");
  const [allowDrag, setAllowDrag] = useState<boolean>(true);
  const [plotLayout, setPlotLayout] = useState<any>({
    width: 220 - 20,
    height: 200 - 60,
  });

  const handleSelectChange = (event: any) => {
    setSelectedOption(event.target.value);
  };

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
          <Form.Select onChange={handleSelectChange}>
            <option value="1">Distance</option>
            <option value="2">plotting/histogram.json</option>
          </Form.Select>
          <Button
            variant="tertiary"
            className=" mx-2 btn btn-outline-secondary"
          >
            <IoDuplicate />
          </Button>
          {/* TODO: tooltip */}
          <Button variant="close" className="mx-2" />
        </Card.Header>
        <Card.Body style={{ padding: 0 }}>
          {selectedOption === "1" && (
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
          )}
          {selectedOption === "2" && (
            <Plot
              data={[
                {
                  x: [1, 2, 3, 4],
                  y: [10, 15, 13, 17],
                  mode: "markers",
                  type: "scatter",
                },
                {
                  x: [2, 3, 4, 5],
                  y: [16, 5, 11, 9],
                  mode: "lines",
                  type: "scatter",
                },
                {
                  x: [1, 2, 3, 4],
                  y: [12, 9, 15, 12],
                  mode: "lines+markers",
                  type: "scatter",
                },
              ]}
              layout={plotLayout}
              onBeforeHover={() => setAllowDrag(false)}
              onUnhover={() => setAllowDrag(true)}
            />
          )}
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
