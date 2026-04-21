import { useFilesystemProviders } from "./useFilesystemProviders";

/**
 * Returns ``true`` when the current room has at least one provider of
 * category ``"filesystem"``.
 */
export function useHasFilesystemProviders(): boolean {
	const { data } = useFilesystemProviders();
	return (data?.length ?? 0) > 0;
}
