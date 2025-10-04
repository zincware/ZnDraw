import * as THREE from 'three';
import { useQueries } from '@tanstack/react-query';
import { getFrameDataOptions } from '../../hooks/useTrajectoryData';
import { useAppStore } from '../../store';
import { useRef, useMemo } from 'react';
import { useFrame } from '@react-three/fiber';
import { BufferGeometryUtils } from "three/examples/jsm/Addons.js";
import { renderMaterial } from './materials';

// Props interface for the dynamic keys
type StaticValue = number | number[] | number[][];
type DataProp = string | StaticValue;

interface ArrowProps {
    start: DataProp;
    direction: DataProp;
    color: DataProp;
    radius: DataProp;
    scale: DataProp;
    material: string;
}

// Reusable THREE objects
const _startVec = new THREE.Vector3();
// No need for _endVec if we use direction directly
const _direction = new THREE.Vector3();
// The default orientation of the arrow geometry is along the +Y axis
const _arrowUp = new THREE.Vector3(0, 1, 0);
const _quaternion = new THREE.Quaternion();
const _scaleVec = new THREE.Vector3();
const _matrix = new THREE.Matrix4();
const _color = new THREE.Color();

function createArrowMesh() {
    const cylinderRadius = 0.04;
    const cylinderHeight = 0.6;
    const coneRadius = 0.1;
    const coneHeight = 0.4;

    const cylinderGeometry = new THREE.CylinderGeometry(
        cylinderRadius,
        cylinderRadius,
        cylinderHeight,
        32,
    );
    const coneGeometry = new THREE.ConeGeometry(coneRadius, coneHeight, 32);

    // Position geometries so the base is at (0,0,0) and it extends up to a total height of 1.0
    cylinderGeometry.translate(0, cylinderHeight / 2, 0);
    coneGeometry.translate(0, cylinderHeight + coneHeight / 2, 0);

    const arrowGeometry = BufferGeometryUtils.mergeGeometries([
        cylinderGeometry,
        coneGeometry,
    ]);

    // This geometry now has a total height of 1.0 and points along the +Y axis.
    return arrowGeometry;
}

export default function Arrow({ start, direction, color, radius, scale, material }: ArrowProps) {
    const instancedMeshRef = useRef<THREE.InstancedMesh | null>(null);
    const { currentFrame, roomId, clientId, selection } = useAppStore();
    const lastGoodFrameData = useRef<any>(null);

    const selectionSet = useMemo(() => {
        return selection ? new Set(selection) : null;
    }, [selection]);

    const { dynamicKeys, staticValues } = useMemo(() => {
        const props = { start, direction, color, radius, scale };
        const dynamicKeys: { [key: string]: string } = {};
        const staticValues: { [key: string]: StaticValue } = {};

        for (const key in props) {
            const value = props[key as keyof ArrowProps];
            if (typeof value === 'string') {
                dynamicKeys[key] = value;
            } else {
                staticValues[key] = value;
            }
        }
        return { dynamicKeys, staticValues };
    }, [start, direction, color, radius, scale]);

    const requiredKeys = Object.values(dynamicKeys);

    const queries = useMemo(() => {
        if (!roomId) {
            return [];
        }
        return requiredKeys.map(key => getFrameDataOptions(roomId, currentFrame, key));
    }, [currentFrame, roomId, requiredKeys]);

    const queryResults = useQueries({ queries });

    // REWRITTEN: The data combination logic
    const { frameData, isFetching } = useMemo(() => {
        const isFetching = queryResults.some(result => result.isFetching || result.isPlaceholderData);

        // 1. Get fetched data
        const fetchedData: { [key: string]: any[] } = {};
        const propToKeyMap = Object.entries(dynamicKeys);
        for (let i = 0; i < requiredKeys.length; i++) {
            const key = requiredKeys[i];
            const result = queryResults[i];
            if (result?.data?.data) {
                // Find which prop this key belongs to (e.g., "arrays.positions" -> "start")
                const propName = propToKeyMap.find(([prop, k]) => k === key)?.[0];
                if (propName) {
                    fetchedData[propName] = result.data.data;
                }
            }
        }

        // 2. Determine the instance count
        let count = 0;
        const firstFetchedProp = Object.values(fetchedData)[0];
        if (firstFetchedProp) {
            // Count from dynamic data (assuming 3 components per instance e.g. [x,y,z,x,y,z...])
            count = firstFetchedProp.length / 3;
        } else if (staticValues.start && Array.isArray(staticValues.start[0])) {
            // Count from static 'start' array
            count = staticValues.start.length;
        } else if (staticValues.direction && Array.isArray(staticValues.direction[0])) {
            // Or count from static 'direction' array
            count = staticValues.direction.length;
        }
        if (count === 0) return { isFetching, frameData: null };

        // 3. Combine fetched and static data into a final structure
        const combinedData: { [key: string]: any } = { count };
        const allPropNames = ['start', 'direction', 'color', 'radius', "scale"];

        for (const prop of allPropNames) {
            if (fetchedData[prop]) {
                // Use fetched data
                combinedData[prop] = fetchedData[prop];
            } else if (staticValues[prop] !== undefined) {
                // Use static data, expanding it if it's a uniform value
                const staticVal = staticValues[prop];
                if (Array.isArray(staticVal) && (Array.isArray(staticVal[0]) || staticVal.length === count * 3)) {
                    // Per-instance array (e.g., [[x,y,z], ...] or flat [x,y,z,x,y,z...])
                    combinedData[prop] = Array.isArray(staticVal[0]) ? staticVal.flat() : staticVal;
                } else {
                    // Uniform value (e.g., 5 or [1,0,0]) - expand it for all instances
                    const components = prop === 'radius' || prop === 'scale' ? 1 : 3;
                    const arr = new Float32Array(count * components);
                    for (let i = 0; i < count; i++) {
                        if (components === 1) arr[i] = staticVal as number;
                        else arr.set(staticVal as number[], i * 3);
                    }
                    combinedData[prop] = arr;
                }
            } else {
                // This prop was not provided, return null
                return { isFetching, frameData: null };
            }
        }
        // Rename 'start' to 'startKey' etc. to match useFrame expectations
        const finalFrameData = {
            ["startKey"]: combinedData.start,
            ["directionKey"]: combinedData.direction,
            ["colorKey"]: combinedData.color,
            ["radiusKey"]: combinedData.radius,
            ["scaleKey"]: combinedData.scale,
            count: combinedData.count
        }

        return { isFetching, frameData: finalFrameData };

    }, [queryResults, staticValues, dynamicKeys, requiredKeys]);

    const dataToRender = frameData || lastGoodFrameData.current;

    const geometry = useMemo(createArrowMesh, []);

    useFrame(() => {
        if (!instancedMeshRef.current || !dataToRender || isFetching) {
            return;
        }

        const mesh = instancedMeshRef.current;

        const positions = dataToRender["startKey"];
        const directions = dataToRender["directionKey"];
        const colors = dataToRender["colorKey"];
        const radii = dataToRender["radiusKey"];
        const scales = dataToRender["scaleKey"];
        const { count } = dataToRender;

        // MODIFICATION: Ensure directions data is also present
        if (!positions || !directions || !colors || !radii || !scales || !count) {
            return;
        }
        console.log("Rendering arrows with scale:", scales);

        for (let i = 0; i < count; i++) {
            const i3 = i * 3;

            // Set the starting position of the arrow
            _startVec.set(positions[i3], positions[i3 + 1], positions[i3 + 2]);

            // Set the direction vector
            _direction.set(directions[i3], directions[i3 + 1], directions[i3 + 2]);

            // --- START OF NEW LOGIC ---

            // 1. Calculate the length of the arrow from the direction vector's magnitude.
            const length = _direction.length() * scales[i];

            // 2. Set the scale. The Y-axis corresponds to the arrow's height (length),
            //    while X and Z correspond to its width (radius).
            const radius = radii[i];
            _scaleVec.set(radius, length, radius);

            // 3. Calculate the rotation quaternion.
            //    This finds the rotation needed to go from our default arrow direction (_arrowUp)
            //    to the new target direction. We normalize the direction vector for this calculation.
            if (length > 0.0001) { // Avoid issues with zero-length vectors
                _quaternion.setFromUnitVectors(_arrowUp, _direction.normalize());
            } else {
                // If length is zero, no rotation is needed, and scale on Y is zero, hiding it.
                _quaternion.identity();
            }

            // 4. Compose the final transformation matrix from position, rotation, and scale.
            //    This is more efficient and robust than setting them individually.
            _matrix.compose(_startVec, _quaternion, _scaleVec);

            // --- END OF NEW LOGIC ---

            mesh.setMatrixAt(i, _matrix);

            if (selectionSet && selectionSet.has(i)) {
                _color.setRGB(1.0, 0.75, 0.8); // Pink color
            } else {
                _color.setRGB(colors[i3], colors[i3 + 1], colors[i3 + 2]);
            }
            mesh.setColorAt(i, _color);
        }

        mesh.instanceMatrix.needsUpdate = true;
        if (mesh.instanceColor) {
            mesh.instanceColor.needsUpdate = true;
        }
    });

    if (!clientId || !roomId || !dataToRender) {
        return null;
    }

    if (frameData) {
        lastGoodFrameData.current = frameData;
    }

    return (
        <instancedMesh
            key={dataToRender.count} // Unique key based on count
            ref={instancedMeshRef}
            args={[geometry, undefined, dataToRender.count]}
        >
            {renderMaterial(material)}
        </instancedMesh>
    );
}