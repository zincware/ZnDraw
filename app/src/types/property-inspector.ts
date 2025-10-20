/**
 * Property Inspector types for per-particle and global property display.
 */

export interface PropertyMetadata {
  dtype: string;
  shape?: number[];
  type: "array" | "scalar";
}

export interface PropertyInfo {
  key: string;
  metadata: PropertyMetadata;
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
