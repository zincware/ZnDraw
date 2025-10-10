import { useRef } from "react";
import { Box } from "@react-three/drei";
import { useFrame } from "@react-three/fiber";
import { useColorScheme, useTheme } from "@mui/material/styles";
import * as THREE from "three";

export default function Crosshair() {
    const groupRef = useRef<THREE.Object3D>(null);
    const shortDimension = 0.05;
    const longDimension = 0.5;
    const { mode } = useColorScheme();
    const theme = useTheme();
    
    const color = mode === "dark" ? theme.palette.primary.light : theme.palette.primary.dark;

    useFrame(({ controls }) => {
        if (groupRef.current) {
            if (controls && 'target' in controls) {
                // @ts-ignore - Drei controls add the target property dynamically
                groupRef.current.position.copy(controls.target);
            } 
        }
    });

    return (
        <group ref={groupRef}>
            {/* X axis box */}
            <Box args={[longDimension, shortDimension, shortDimension]}>
                <meshStandardMaterial
                    color={color}
                />
            </Box>

            {/* Y axis box */}
            <Box args={[shortDimension, longDimension, shortDimension]}>
                <meshStandardMaterial
                    color={color}
                />
            </Box>

            {/* Z axis box */}
            <Box args={[shortDimension, shortDimension, longDimension]}>
                <meshStandardMaterial
                    color={color}
                />
            </Box>
        </group>
    );
}