import * as THREE from "three";
import { useQueries } from "@tanstack/react-query";
import { getFrameDataOptions } from "../../hooks/useTrajectoryData";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState } from "react";
import { useFrame } from "@react-three/fiber";
import { renderMaterial } from "./materials";

// Interface for interaction settings (hovering, selecting)
interface InteractionSettings {
    enabled: boolean
    color: string | null
    opacity: number
}

// Interface for the component's data prop
// Position, color, and radius can be a string key (for fetching) or a direct array/value.
interface SphereData {
    position: string | number[][] | number[];
    color: string | number[][] | number[];
    radius: string | number[] | number;
    material: string;
    resolution: number;
    scale: number;
    opacity: number;
    selecting: InteractionSettings;
    hovering: InteractionSettings;
}

// Reusable THREE objects to avoid creating them in the render loop, improving performance.
const positionVec = new THREE.Vector3();
const scaleVec = new THREE.Vector3();
const matrix = new THREE.Matrix4();
const tempColor = new THREE.Color();

// Constants for interaction scaling effects
const HOVER_SCALE = 1.25;
const SELECTION_SCALE = 1.01;

export default function Sphere({ data }: { data: SphereData }) {
    // Destructure props, renaming to avoid naming conflicts within the component
    const {
        position: positionProp,
        color: colorProp,
        radius: radiusProp,
        material,
        resolution,
        scale,
        selecting,
        hovering,
    } = data;

    // Refs for the THREE.js mesh objects
    const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
    const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
    const hoverMeshRef = useRef<THREE.Mesh | null>(null);

    // State for tracking the currently hovered particle instance
    const [hoveredInstanceId, setHoveredInstanceId] = useState<number | null>(null);

    // Zustand store for global state like current frame, selection, etc.
    const { currentFrame, roomId, clientId, selection, updateSelection } = useAppStore();
    const lastGoodFrameData = useRef<any>(null);

    // Memoize the set of selected indices for efficient lookups
    const selectionSet = useMemo(() => new Set(selection || []), [selection]);
    const selectedIndices: number[] = useMemo(() => Array.from(selectionSet), [selectionSet]);

    // Particle geometry and scaling settings
    const particleResolution = resolution || 8;
    const particleScale = scale || 1.0;

    // Step 1: Determine which properties are string keys that require fetching.
    const keysToFetch = useMemo(() => {
        const keys: string[] = [];
        if (typeof positionProp === 'string') keys.push(positionProp);
        if (typeof colorProp === 'string') keys.push(colorProp);
        if (typeof radiusProp === 'string') keys.push(radiusProp);
        return keys;
    }, [positionProp, colorProp, radiusProp]);

    // Step 2: Set up TanStack queries for the keys that need to be fetched.
    const queries = useMemo(() => {
        if (!roomId) return [];
        return keysToFetch.map((key) =>
            getFrameDataOptions(roomId, currentFrame, key),
        );
    }, [currentFrame, roomId, keysToFetch]);

    const queryResults = useQueries({ queries });

    // Step 3: Process and normalize local and fetched data into a consistent format.
    const { processedData, isFetching } = useMemo(() => {
        const isQueryFetching = queryResults.some(
            (result) => result.isFetching || result.isPlaceholderData,
        );

        // If we are fetching data, ensure all queries have succeeded before proceeding.
        const allQueriesSucceeded = queryResults.every((result) => result.isSuccess);
        if (keysToFetch.length > 0 && !allQueriesSucceeded) {
            return { isFetching: isQueryFetching, processedData: null };
        }

        // Create a map of fetched data for easy access.
        const fetchedDataMap = new Map(
            queryResults.map((result, i) => [keysToFetch[i], result.data?.data])
        );

        // --- Normalize Position and determine particle count ---
        let count = 0;
        let finalPositions: number[] = [];
        if (typeof positionProp === 'string') {
            const remoteData = fetchedDataMap.get(positionProp) || [];
            count = remoteData.length / 3;
            finalPositions = remoteData;
        } else {
            if (positionProp.length > 0 && Array.isArray(positionProp[0])) {
                // Data is number[][]: [[x,y,z], ...]
                count = positionProp.length;
                finalPositions = (positionProp as number[][]).flat();
            } else if (positionProp.length > 0) {
                // Data is number[]: [x,y,z]
                count = 1;
                finalPositions = positionProp as number[];
            }
        }

        if (count === 0) {
            return { isFetching: isQueryFetching, processedData: null };
        }

        // --- Normalize Color ---
        let finalColors: number[];
        if (typeof colorProp === 'string') {
            finalColors = fetchedDataMap.get(colorProp) || [];
        } else {
            if (colorProp.length > 0 && Array.isArray(colorProp[0])) {
                // Data is number[][]: [[r,g,b], ...]
                finalColors = (colorProp as number[][]).flat();
            } else {
                // Data is number[]: [r,g,b], repeat for all particles
                finalColors = new Array(count).fill(colorProp).flat();
            }
        }

        // --- Normalize Radius ---
        let finalRadii: number[];
        if (typeof radiusProp === 'string') {
            finalRadii = fetchedDataMap.get(radiusProp) || [];
        } else {
            if (Array.isArray(radiusProp)) {
                // Data is number[]: [r1, r2, ...]
                finalRadii = radiusProp as number[];
            } else {
                // Data is a single number, repeat for all particles
                finalRadii = new Array(count).fill(radiusProp);
            }
        }

        // Final sanity check to ensure data arrays are consistent
        if (finalPositions.length / 3 !== count || finalColors.length / 3 !== count || finalRadii.length !== count) {
            console.error("Sphere data arrays have inconsistent lengths.");
            return { isFetching: isQueryFetching, processedData: null };
        }
        console.log("Sphere processed data count:", { count, finalPositions, finalColors, finalRadii });

        return {
            isFetching: isQueryFetching,
            processedData: { positions: finalPositions, colors: finalColors, radii: finalRadii, count },
        };
    }, [queryResults, keysToFetch, positionProp, colorProp, radiusProp]);

    // Use the newly processed data, or fallback to the last known good data during fetches.
    const dataToRender = processedData || lastGoodFrameData.current;
    if (processedData) {
        lastGoodFrameData.current = processedData;
    }

    // useFrame hook runs on every rendered frame, ideal for updating our THREE meshes.
    useFrame(() => {
        if (!mainMeshRef.current || !dataToRender || isFetching) {
            return;
        }

        const mainMesh = mainMeshRef.current;
        const selectionMesh = selectionMeshRef.current;
        const hoverMesh = hoverMeshRef.current;
        const { positions, colors, radii, count } = dataToRender;

        // Update the main instanced mesh for all particles
        for (let i = 0; i < count; i++) {
            const i3 = i * 3;
            positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
            const r = radii[i] * particleScale;
            scaleVec.set(r, r, r);
            matrix.identity().setPosition(positionVec).scale(scaleVec);
            mainMesh.setMatrixAt(i, matrix);
            tempColor.setRGB(colors[i3], colors[i3 + 1], colors[i3 + 2]);
            mainMesh.setColorAt(i, tempColor);
        }

        if (selectionMesh) {
            selectedIndices.forEach((id, index) => {
                if (id >= count) return;
                const i3 = id * 3;
                positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
                const r = radii[id] * particleScale * SELECTION_SCALE;
                scaleVec.set(r, r, r);
                matrix.identity().setPosition(positionVec).scale(scaleVec);
                selectionMesh.setMatrixAt(index, matrix);
            });
        }

        if (hoverMesh) {
            if (hovering.enabled && hoveredInstanceId !== null && hoveredInstanceId < count) {
                hoverMesh.visible = true;
                const i3 = hoveredInstanceId * 3;
                positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
                const r = radii[hoveredInstanceId] * particleScale * HOVER_SCALE;
                hoverMesh.position.copy(positionVec);
                hoverMesh.scale.set(r, r, r);
            } else {
                hoverMesh.visible = false;
            }
        }

        // Set flags to tell THREE.js to update the buffers.
        mainMesh.count = count;
        mainMesh.instanceMatrix.needsUpdate = true;
        if (mainMesh.instanceColor) mainMesh.instanceColor.needsUpdate = true;

        // Conditionally flag updates for selection mesh
        if (selectionMesh) {
            selectionMesh.count = selectedIndices.length;
            selectionMesh.instanceMatrix.needsUpdate = true;
        }
    });

    if (!clientId || !roomId || !dataToRender) {
        return null;
    }

    // --- Event Handlers ---
    const onClickHandler = (event: any) => {
        if (event.detail !== 1 || event.instanceId === undefined) return;
        event.stopPropagation();
        updateSelection(event.instanceId, event.shiftKey);
    };

    const onPointerMoveHandler = (event: any) => {
        if (event.instanceId === undefined) return;
        event.stopPropagation();
        setHoveredInstanceId(event.instanceId);
    };

    const onPointerOutHandler = () => setHoveredInstanceId(null);

    return (
        <group>
            {/* Main mesh for all particles */}
            <instancedMesh
                key={dataToRender.count}
                ref={mainMeshRef}
                args={[undefined, undefined, dataToRender.count]}
                onClick={selecting.enabled ? onClickHandler : undefined}
                onPointerMove={hovering.enabled ? onPointerMoveHandler : undefined}
                onPointerOut={hovering.enabled ? onPointerOutHandler : undefined}
            >
                <sphereGeometry args={[1, particleResolution, particleResolution]} />
                {renderMaterial(material, data.opacity)}
            </instancedMesh>

            {/* Instanced mesh for selected particles */}
            {selecting.enabled && (
                <instancedMesh
                    key={`selection-${selectedIndices.length}`}
                    ref={selectionMeshRef}
                    args={[undefined, undefined, selectedIndices.length]}
                >
                    <sphereGeometry args={[1, particleResolution, particleResolution]} />
                    <meshBasicMaterial
                        side={THREE.FrontSide}
                        transparent
                        opacity={selecting.opacity}
                        color={selecting.color || "#FFFF00"}
                    />
                </instancedMesh>
            )}

            {/* Single mesh for the hovered particle */}
            {hovering.enabled && (
                <mesh ref={hoverMeshRef} visible={false}>
                    <sphereGeometry args={[1, particleResolution, particleResolution]} />
                    <meshBasicMaterial
                        side={THREE.BackSide}
                        transparent
                        opacity={hovering.opacity}
                        color={hovering.color || "#00FFFF"}
                    />
                </mesh>
            )}
        </group>
    );
}
