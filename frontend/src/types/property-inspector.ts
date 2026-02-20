/**
 * Property Inspector types for per-particle and global property display.
 */
import type { PropertyMeta } from "../myapi/client";

export interface PropertyInfo {
	key: string;
	metadata: PropertyMeta;
	category: "per-particle" | "global";
	prefix: string; // "arrays", "calc", "info", etc.
	enabled: boolean;
}

export interface PropertyValue {
	key: string;
	value: any; // Actual data from API
	formattedValue: string; // For display
	isLoading: boolean;
	error?: string;
}

export interface PropertyInspectorState {
	selectedProperties: Set<string>;
	expandedCategories: Set<string>;
	searchQuery: string;
	sortBy: "name" | "category" | "type";
	filterBy: "all" | "per-particle" | "global";
}
