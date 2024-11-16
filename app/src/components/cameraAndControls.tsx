import { useCallback, useEffect, useMemo, useRef } from "react";
import {
  PerspectiveCamera,
  OrthographicCamera,
  TrackballControls,
  OrbitControls,
} from "@react-three/drei";
import * as THREE from "three";
import { debounce } from "lodash";
import { useCentroid } from "./particlesEditor";

type CameraAndControls = {
  camera: THREE.Vector3;
  target: THREE.Vector3;
};

type CameraAndControlsProps = {
    roomConfig: any;
    cameraAndControls: CameraAndControls;
    setCameraAndControls: React.Dispatch<React.SetStateAction<CameraAndControls>>;
    currentFrame: any;
    selectedIds: Set<number>;
  };

  const CameraAndControls: React.FC<CameraAndControlsProps> = ({
    roomConfig,
    cameraAndControls,
    setCameraAndControls,
    currentFrame,
    selectedIds
  }) => {
    const cameraRef = useRef<any>(null);
    const controlsRef = useRef<any>(null);
    const centroid = useCentroid({ frame: currentFrame, selectedIds: selectedIds });

    const controlsOnEndFn = useCallback(
        debounce(() => {
            if (!cameraRef.current || !controlsRef.current) {
                return;
            }
            setCameraAndControls({
                camera: cameraRef.current.position,
                target: controlsRef.current.target,
            });
        }, 100),
        []
    );

    // if the camera positions and target positions is default, adapt them to the scene
    useEffect(() => {
        if (cameraAndControls.camera.equals(new THREE.Vector3(0, 0, 0)) && cameraAndControls.target.equals(new THREE.Vector3(0, 0, 0))) {
            if (!centroid) {
                return;
            }
            if (currentFrame.positions.length === 0) {
                return;
            }
            if (cameraRef.current === null) {
                return;
            }
            if (centroid.equals(new THREE.Vector3(0, 0, 0))) {
                return;
            }
            // check that centroid is not NaN
            if (isNaN(centroid.x) || isNaN(centroid.y) || isNaN(centroid.z)) {
                return;
            }
            // now calculate the camera positions 

            // Compute the bounding sphere radius
            let maxDistance = 0;
            currentFrame.positions.forEach((x) => {
                maxDistance = Math.max(maxDistance, x.distanceTo(centroid));
            });

            const fov = (cameraRef.current.fov * Math.PI) / 180; // Convert FOV to radians
            const distance = maxDistance / Math.tan(fov / 2);

            console.log("cameraAndControls default", cameraAndControls);
            console.log("centroid", centroid);
            setCameraAndControls({
                camera: new THREE.Vector3(distance, distance, distance),
                target: centroid
            });
        }
    }, [centroid, currentFrame.positions]);

    useEffect(() => {
        if (!cameraRef.current || !controlsRef.current ) {
            return;
        }
        console.log("cameraAndControls update", cameraAndControls);
        // controlsRef.current.enabled = false;
        cameraRef.current.position.copy(cameraAndControls.camera);
        controlsRef.current.target.copy(cameraAndControls.target);
        controlsRef.current.update();
        // controlsRef.current.enabled = true;
    }, [cameraAndControls]);

  return (
    <>  
        {roomConfig.scene.camera === "OrthographicCamera" && (
            <OrthographicCamera makeDefault ref={cameraRef}>
                <pointLight
                    decay={0}
                    intensity={Math.PI / 2}
                />
            </OrthographicCamera>
        )}
        {roomConfig.scene.camera === "PerspectiveCamera" && (
            <PerspectiveCamera makeDefault ref={cameraRef}>
                <pointLight
                    decay={0}
                    intensity={Math.PI / 2}
                />
            </PerspectiveCamera>
        )}
        <OrbitControls makeDefault onEnd={controlsOnEndFn} enableDamping={false} ref={controlsRef}/>            
    </>
  );
};

export default CameraAndControls;