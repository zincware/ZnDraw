import { useEffect, useState } from "react";

import * as THREE from "three";

export const useColorMode = (): [string, () => void] => {
	const [colorMode, setColorMode] = useState<string>("light");

	useEffect(() => {
		const theme =
			localStorage.getItem("theme") ||
			(window.matchMedia("(prefers-color-scheme: dark)").matches
				? "dark"
				: "light");
		setTheme(theme, setColorMode);
	}, []);

	const handleColorMode = () => {
		const newColorMode = colorMode === "light" ? "dark" : "light";
		setTheme(newColorMode, setColorMode);
		localStorage.setItem("theme", newColorMode);
	};

	return [colorMode, handleColorMode];
};

const setTheme = (
	theme: string,
	setColorMode: (mode: string) => void,
): void => {
	document.documentElement.setAttribute("data-bs-theme", theme);
	setColorMode(theme);
};

export type HSLColor = [number, number, number];
export type ColorRange = [number, number];

export const interpolateColor = (
	colors: HSLColor[],
	range: ColorRange,
	value: number,
): THREE.Color => {
	const [min, max] = range;

	// Clamp the value to the range
	if (value <= min) {
		return new THREE.Color().setHSL(...colors[0]);
	}
	if (value >= max) {
		return new THREE.Color().setHSL(...colors[colors.length - 1]);
	}

	// Normalize the value within the range
	const normalizedValue = (value - min) / (max - min);

	// Calculate the exact position within the colors array
	const scaledValue = normalizedValue * (colors.length - 1);
	const lowerIndex = Math.floor(scaledValue);
	const upperIndex = Math.ceil(scaledValue);

	const lowerColor = colors[lowerIndex];
	const upperColor = colors[upperIndex];
	const t = scaledValue - lowerIndex;

	// Interpolate between lowerColor and upperColor
	const h = THREE.MathUtils.lerp(lowerColor[0], upperColor[0], t);
	const s = THREE.MathUtils.lerp(lowerColor[1], upperColor[1], t);
	const l = THREE.MathUtils.lerp(lowerColor[2], upperColor[2], t);

	return new THREE.Color().setHSL(h, s, l);
};

// Define the type for your state
export type IndicesState = {
	active: boolean;
	indices: Set<number>; // or Set<string> if you have string indices
};
