import * as THREE from "three";
import { useQueries } from "@tanstack/react-query";
import { getFrameDataOptions } from "../../hooks/useTrajectoryData";
import { useAppStore } from "../../store";
import { useRef, useMemo, useState } from "react";
import { useFrame } from "@react-three/fiber";
import { renderMaterial } from "./materials";

interface InteractionSettings {
    enabled: boolean
    color: string | null
    opacity: number
}

interface SphereData {
    position: string;
    color: string;
    radius: string;
    material: string;
    resolution: number;
    scale: number;
    selecting: InteractionSettings;
    hovering: InteractionSettings;
}

// Reusable vectors and objects to avoid creating them in the loop
const positionVec = new THREE.Vector3();
const scaleVec = new THREE.Vector3();
const matrix = new THREE.Matrix4();
const color = new THREE.Color();

const HOVER_SCALE = 1.25;
const SELECTION_SCALE = 1.01;

export default function Sphere({ data }: { data: SphereData }) {
    const {
        position: positionKey,
        color: colorKey,
        radius: radiusKey,
        material,
        resolution,
        scale,
        selecting,
        hovering,
    } = data;

    const mainMeshRef = useRef<THREE.InstancedMesh | null>(null);
    const selectionMeshRef = useRef<THREE.InstancedMesh | null>(null);
    const hoverMeshRef = useRef<THREE.Mesh | null>(null);

    const [hoveredInstanceId, setHoveredInstanceId] = useState<number | null>(null);

    const { currentFrame, roomId, clientId, selection, updateSelection } = useAppStore();
    const lastGoodFrameData = useRef<any>(null);

    const selectionSet = useMemo(() => {
        return selection ? new Set(selection) : new Set();
    }, [selection]);

    const selectedIndices: number[] = useMemo(() => Array.from(selectionSet), [selectionSet]);

    const particleResolution = resolution || 8;
    const particleScale = scale || 1.0;

    const requiredKeys = useMemo(
        () => [positionKey, colorKey, radiusKey],
        [positionKey, colorKey, radiusKey],
    );

    const queries = useMemo(() => {
        if (!roomId) {
            return [];
        }
        return requiredKeys.map((key) =>
            getFrameDataOptions(roomId, currentFrame, key),
        );
    }, [currentFrame, roomId, requiredKeys]);

    const queryResults = useQueries({ queries });

    const { frameData, isFetching } = useMemo(() => {
        const isFetching = queryResults.some(
            (result) => result.isFetching || result.isPlaceholderData,
        );
        const allDataPresent = queryResults.every((result) => result.data);

        if (!allDataPresent || queryResults.length === 0) {
            return { isFetching, frameData: null };
        }

        const firstSuccessfulResult = queryResults.find(
            (result) => result.isSuccess,
        );
        const combinedData: { [key: string]: any } = {};
        let isComplete = true;

        for (let i = 0; i < requiredKeys.length; i++) {
            const key = requiredKeys[i];
            const result = queryResults[i];
            if (result?.data?.data) {
                combinedData[key] = result.data.data;
            } else {
                isComplete = false;
                break;
            }
        }

        if (!isComplete) {
            return { isFetching, frameData: null };
        }

        combinedData.count = firstSuccessfulResult?.data?.shape[0] || 0;
        return { isFetching, frameData: combinedData };
    }, [queryResults, requiredKeys]);

    const dataToRender = frameData || lastGoodFrameData.current;

    useFrame(() => {
        if (!mainMeshRef.current || !selectionMeshRef.current || !hoverMeshRef.current || !dataToRender || isFetching) {
            return;
        }

        const mainMesh = mainMeshRef.current;
        const selectionMesh = selectionMeshRef.current;
        const hoverMesh = hoverMeshRef.current;

        const positions = dataToRender[positionKey];
        const colors = dataToRender[colorKey];
        const radii = dataToRender[radiusKey];
        const { count } = dataToRender;

        if (!positions || !colors || !radii || !count) {
            return;
        }

        // 1. Update the main instanced mesh
        for (let i = 0; i < count; i++) {
            const i3 = i * 3;
            positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);

              const r = radii[i] * particleScale;
              scaleVec.set(r, r, r);

            matrix.identity().setPosition(positionVec).scale(scaleVec);
            mainMesh.setMatrixAt(i, matrix);
            color.setRGB(colors[i3], colors[i3 + 1], colors[i3 + 2]);
            mainMesh.setColorAt(i, color);
        }

        // 2. Update the selection mesh
        selectedIndices.forEach((id, index) => {
            const i3 = id * 3;
            positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
            const r = radii[id] * particleScale * SELECTION_SCALE;
            scaleVec.set(r, r, r);

            matrix.identity().setPosition(positionVec).scale(scaleVec);
            selectionMesh.setMatrixAt(index, matrix);
        });

        // 3. Update the hover sphere
        if (hovering.enabled && hoveredInstanceId !== null && positions[hoveredInstanceId * 3] !== undefined) {
            hoverMesh.visible = true;
            const i3 = hoveredInstanceId * 3;
            positionVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);
            const r = radii[hoveredInstanceId] * particleScale * HOVER_SCALE;

            hoverMesh.position.copy(positionVec);
            hoverMesh.scale.set(r, r, r);
        } else {
            hoverMesh.visible = false;
        }

        mainMesh.instanceMatrix.needsUpdate = true;
        if (mainMesh.instanceColor) {
            mainMesh.instanceColor.needsUpdate = true;
        }
        selectionMesh.instanceMatrix.needsUpdate = true;
    });

    if (!clientId || !roomId || !dataToRender) {
        return null;
    }

    if (frameData) {
        lastGoodFrameData.current = frameData;
    }

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

    const onPointerOutHandler = () => {
        setHoveredInstanceId(null);
    };

    return (
        <group>
            {/* Main mesh for all non-selected particles */}
            <instancedMesh
                key={dataToRender.count}
                ref={mainMeshRef}
                args={[undefined, undefined, dataToRender.count]}
                onClick={selecting.enabled ? onClickHandler : undefined}
                onPointerMove={hovering.enabled ? onPointerMoveHandler : undefined}
                onPointerOut={hovering.enabled ? onPointerOutHandler : undefined}
            >
                <sphereGeometry args={[1, particleResolution, particleResolution]} />
                {renderMaterial(material)}
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