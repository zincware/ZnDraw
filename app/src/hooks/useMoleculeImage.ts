// src/hooks/useMoleculeImage.ts
// Custom hook for fetching and managing molecule images from SMILES notation
// Uses react-query for efficient caching and deduplication

import { useState, useEffect } from "react";
import { useQuery } from "@tanstack/react-query";
import { useColorScheme } from "@mui/material/styles";
import { convertMoleculeToImage } from "../myapi/client";

interface MoleculeImageState {
	image: string | null;
	loading: boolean;
	error: string | null;
}

/**
 * Custom hook for fetching molecule images from SMILES notation.
 * Uses react-query for automatic caching, deduplication, and stale data management.
 * Automatically uses dark mode based on the current theme.
 *
 * @param smiles - SMILES notation string
 * @param throttleMs - Optional throttle delay in milliseconds (default: 0)
 * @returns Object containing image URL, loading state, and error message
 */
export const useMoleculeImage = (
	smiles: string,
	throttleMs: number = 0,
): MoleculeImageState => {
	const [debouncedSmiles, setDebouncedSmiles] = useState(smiles);
	const { mode, systemMode } = useColorScheme();
	// Only use dark mode when explicitly "dark" or when system resolves to "dark"
	// Default to light mode for undefined/initial state to avoid grey bonds
	const isDarkMode =
		mode === "dark" || (mode === "system" && systemMode === "dark");

	// Debounce/throttle the SMILES input
	useEffect(() => {
		if (throttleMs > 0) {
			const timeoutId = setTimeout(() => {
				setDebouncedSmiles(smiles);
			}, throttleMs);

			return () => clearTimeout(timeoutId);
		} else {
			setDebouncedSmiles(smiles);
		}
	}, [smiles, throttleMs]);

	// Use react-query for fetching with automatic caching
	const {
		data: image,
		isLoading: loading,
		error: queryError,
	} = useQuery({
		queryKey: ["molecule-image", debouncedSmiles, isDarkMode],
		queryFn: async () => {
			if (!debouncedSmiles || debouncedSmiles.trim() === "") {
				return null;
			}

			const response = await convertMoleculeToImage({
				type: "smiles",
				data: debouncedSmiles,
				dark: isDarkMode,
			});

			return response.image;
		},
		enabled: Boolean(debouncedSmiles && debouncedSmiles.trim() !== ""),
		staleTime: Infinity, // Molecule images don't change, cache forever
		gcTime: 1000 * 60 * 30, // Keep in cache for 30 minutes
		retry: 1, // Only retry once on failure
	});

	// Format error message
	const error = queryError
		? (queryError as any)?.response?.data?.error || "Failed to generate image"
		: null;

	return {
		image: image || null,
		loading,
		error,
	};
};
