import type React from "react";
import {
    type ReactNode,
    createContext,
    useEffect,
    useState,
    useContext,
    useMemo,
} from "react";
import { useAppContext } from "./AppContext";
import { useFrameConnection } from "../components/api";
import { Dict } from "znsocket";
import { set } from "lodash";

interface SlowFrameState {
    isSlowFrame: boolean;
    atomsInfo: Record<string, any>;
    atomsArrays: Record<string, any>;
    threshold: number; // Optional threshold for slow frame detection
}

const SlowFrameContext = createContext<SlowFrameState | undefined>(undefined);

export const useSlowFrame = () => {
    const context = useContext(SlowFrameContext);
    if (context === undefined) {
        throw new Error("useSlowFrame must be used within a SlowFrameProvider");
    }
    return context;
};

interface SlowFrameProviderProps {
    children: ReactNode;
    threshold?: number;
}

export const SlowFrameProvider: React.FC<SlowFrameProviderProps> = ({
    children,
    threshold = 150,
}) => {
    const { step, token } = useAppContext();
    const { framesCon } = useFrameConnection(token);

    const [isSlowFrame, setIsSlowFrame] = useState(false);
    const [atomsInfo, setAtomsInfo] = useState<Record<string, any>>({});
    const [atomsArrays, setAtomsArrays] = useState<Record<string, any>>({});

    useEffect(() => {
        const fetchData = async () => {
            try {
                const frameData = await framesCon?.get(step);
                if (!frameData) return;

                // Fetch info and calc concurrently
                const [info, calc, arrays] = await Promise.all([
                    frameData.get("info") as Promise<Dict | null>,
                    frameData.get("calc") as Promise<Dict | null>,
                    frameData.get("arrays") as Promise<Dict | null>,
                ]);

                let newEntries = {};

                if (info) {
                    const infoEntries = await info.entries();
                    newEntries = { ...newEntries, ...Object.fromEntries(infoEntries) };
                }
                if (calc) {
                    const calcEntries = await calc.entries();
                    newEntries = { ...newEntries, ...Object.fromEntries(calcEntries) };
                }
                if (arrays) {
                    const arraysEntries = await arrays.entries();
                    setAtomsArrays(() => ({...Object.fromEntries(arraysEntries) }));
                }
                // iterate over the entries, check if any of them are Dict and if so, convert them to a plain object
                for (const [key, value] of Object.entries(newEntries)) {
                    if (value instanceof Dict) {
                        let plainObject = await value.entries();
                        newEntries[key] = Object.fromEntries(plainObject);
                    }
                }

                setAtomsInfo((prev) => ({ ...prev, ...newEntries }));
                setIsSlowFrame(true);

            } catch (error) {
                console.error("Failed to fetch slow frame data:", error);
            }
        };

        const timer = setTimeout(fetchData, threshold);

        return () => {
            clearTimeout(timer);
            setIsSlowFrame(false);

            setAtomsInfo((prev) => {
                if (Object.keys(prev).length === 0) {
                    return prev; // Return current state to prevent re-render
                }
                return {}; // Otherwise, reset to a new empty object
            });
            setAtomsArrays((prev) => {
                if (Object.keys(prev).length === 0) {
                    return prev; // Return current state to prevent re-render
                }
                return {}; // Otherwise, reset to a new empty object
            });
            // TODO: need to disconnect the framesCon?
        };
    }, [step, threshold, framesCon]);

    useEffect(() => {
        console.log("atomsInfo changed:", atomsInfo);
    }, [atomsInfo]);

    const value = useMemo(() => ({ isSlowFrame, atomsInfo, threshold, atomsArrays }), [isSlowFrame, atomsInfo, threshold, atomsArrays]);

    return (
        <SlowFrameContext.Provider value={value}>
            {children}
        </SlowFrameContext.Provider>
    );
};
