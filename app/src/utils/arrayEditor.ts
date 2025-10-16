/**
 * Utility functions and types for array editing in geometry forms
 */

/**
 * Field type determines how arrays are structured and validated
 */
export type ArrayFieldType =
  | 'position'      // [x, y, z] per instance
  | 'direction'     // [x, y, z] per instance
  | 'color'         // hex color string per instance (e.g., "#FF0000")
  | 'radius'        // single value or per instance
  | 'scale'         // single value or per instance
  | 'rotation'      // [x, y, z] per instance
  | 'size_2d'       // [width, height] per instance (Plane)
  | 'size_3d'       // [width, height, depth] per instance (Box)
  | 'generic';      // unknown type

/**
 * Configuration for different field types
 */
export interface FieldTypeConfig {
  /** Number of columns for this field type (e.g., 3 for [x,y,z]) */
  dimensions: number;
  /** Column labels */
  columnLabels: string[];
  /** Value range for validation [min, max] */
  valueRange?: [number, number];
  /** Whether this field supports single-value mode (applies to all instances) */
  supportsSingleValue: boolean;
  /** Default value for new rows/columns */
  defaultValue: number | string;
  /** Whether this field uses string values (e.g., hex colors) instead of numbers */
  isStringType?: boolean;
}

/**
 * Get configuration for a field type
 */
export function getFieldTypeConfig(fieldType: ArrayFieldType): FieldTypeConfig {
  const configs: Record<ArrayFieldType, FieldTypeConfig> = {
    position: {
      dimensions: 3,
      columnLabels: ['X', 'Y', 'Z'],
      supportsSingleValue: false,
      defaultValue: 0,
    },
    direction: {
      dimensions: 3,
      columnLabels: ['X', 'Y', 'Z'],
      supportsSingleValue: false,
      defaultValue: 0,
    },
    color: {
      dimensions: 1,
      columnLabels: ['Color'],
      supportsSingleValue: true,
      defaultValue: '#FF0000',
      isStringType: true,
    },
    radius: {
      dimensions: 1,
      columnLabels: ['Radius'],
      valueRange: [0, Infinity],
      supportsSingleValue: true,
      defaultValue: 1.0,
    },
    scale: {
      dimensions: 1,
      columnLabels: ['Scale'],
      valueRange: [0, Infinity],
      supportsSingleValue: true,
      defaultValue: 1.0,
    },
    rotation: {
      dimensions: 3,
      columnLabels: ['X', 'Y', 'Z'],
      supportsSingleValue: true,  // Rotation can have a shared value
      defaultValue: 0,
    },
    size_2d: {
      dimensions: 2,
      columnLabels: ['Width', 'Height'],
      valueRange: [0, Infinity],
      supportsSingleValue: true,
      defaultValue: 1.0,
    },
    size_3d: {
      dimensions: 3,
      columnLabels: ['Width', 'Height', 'Depth'],
      valueRange: [0, Infinity],
      supportsSingleValue: true,
      defaultValue: 1.0,
    },
    generic: {
      dimensions: 1,
      columnLabels: ['Value'],
      supportsSingleValue: true,
      defaultValue: 0,
    },
  };

  return configs[fieldType];
}

/**
 * Infer field type from schema property name or path
 * @param path - The field path (e.g., "size", "position", "rotation")
 * @param schema - Optional schema to help determine exact type (e.g., 2D vs 3D size)
 */
export function inferFieldType(path: string, schema?: any): ArrayFieldType {
  const pathLower = path.toLowerCase();

  if (pathLower.includes('position')) return 'position';
  if (pathLower.includes('direction')) return 'direction';
  if (pathLower.includes('color') || pathLower.includes('colour')) return 'color';
  if (pathLower.includes('radius')) return 'radius';
  if (pathLower.includes('scale')) return 'scale';
  if (pathLower.includes('rotation') || pathLower.includes('rotate')) return 'rotation';

  // Special handling for "size" - determine dimensions from schema or default to 3D
  if (pathLower.includes('size')) {
    // Check schema description for hints about dimensions
    if (schema?.description) {
      const desc = schema.description.toLowerCase();
      if (desc.includes('width') && desc.includes('height') && desc.includes('depth')) {
        return 'size_3d';
      }
      if (desc.includes('width') && desc.includes('height') && !desc.includes('depth')) {
        return 'size_2d';
      }
    }
    // Default to 3D (Box is more common)
    return 'size_3d';
  }

  return 'generic';
}

/**
 * Normalize value to array format for editing
 * Handles: number, number[], number[][], string, string[]
 */
export function normalizeToArray(
  value: string | number | (string | number)[] | (string | number)[][],
  fieldType: ArrayFieldType
): (string | number)[][] {
  const config = getFieldTypeConfig(fieldType);

  // For string types (like hex colors)
  if (config.isStringType) {
    // If it's already a 2D array of strings, return as-is
    if (Array.isArray(value) && Array.isArray(value[0])) {
      return value as string[][];
    }

    // If it's a 1D array of strings (e.g., ["#FF0000", "#00FF00"])
    if (Array.isArray(value)) {
      return (value as string[]).map(item => [item]);
    }

    // If it's a single string (e.g., "#FF0000")
    if (typeof value === 'string') {
      return [[value]];
    }

    // Default: return single row with default value
    return [[config.defaultValue as string]];
  }

  // For numeric types
  // If it's already a 2D array, return as-is
  if (Array.isArray(value) && Array.isArray(value[0])) {
    return value as number[][];
  }

  // If it's a 1D array
  if (Array.isArray(value)) {
    const arr = value as number[];
    // If dimensions match, wrap in array (single row)
    if (arr.length === config.dimensions) {
      return [arr];
    }
    // Otherwise, chunk by dimensions
    const result: number[][] = [];
    for (let i = 0; i < arr.length; i += config.dimensions) {
      result.push(arr.slice(i, i + config.dimensions));
    }
    return result;
  }

  // If it's a single number, create single row
  if (typeof value === 'number') {
    return [Array(config.dimensions).fill(value)];
  }

  // Default: return single row with default values
  return [Array(config.dimensions).fill(config.defaultValue as number)];
}

/**
 * Convert 2D array back to appropriate format based on field type and row count
 */
export function denormalizeFromArray(
  arrayValue: (string | number)[][],
  fieldType: ArrayFieldType
): string | number | (string | number)[] | (string | number)[][] {
  const config = getFieldTypeConfig(fieldType);

  // Empty array
  if (arrayValue.length === 0) {
    return config.supportsSingleValue ? config.defaultValue : [];
  }

  // For string types (like hex colors)
  if (config.isStringType) {
    // Single row for single-value-supporting fields
    if (arrayValue.length === 1 && config.supportsSingleValue) {
      return arrayValue[0][0] as string;
    }
    // Multiple rows: return 1D array of strings
    return arrayValue.map(row => row[0] as string);
  }

  // For numeric types
  // Single row for single-value-supporting fields
  if (arrayValue.length === 1 && config.supportsSingleValue) {
    // If it's a 1D field (like radius), return single number
    if (config.dimensions === 1) {
      return arrayValue[0][0] as number;
    }
    // Otherwise return the row array
    return arrayValue[0] as number[];
  }

  // For 1D fields with multiple rows, flatten to 1D array [1, 1, 1]
  if (config.dimensions === 1) {
    return arrayValue.map(row => row[0] as number);
  }

  // Position and direction are ALWAYS per-instance (never single shared value)
  // They must always return 2D array format
  if (fieldType === 'position' || fieldType === 'direction') {
    return arrayValue as number[][];
  }

  // For other multi-dimensional fields (rotation, size):
  // - Single row should return 1D array [x, y, z] (shared value for all instances)
  // - Multiple rows return 2D array [[x,y,z], [x,y,z]] (per-instance values)
  if (arrayValue.length === 1) {
    return arrayValue[0] as number[];
  }

  return arrayValue as number[][];
}

/**
 * Validate array data
 */
export function validateArrayData(
  arrayValue: (string | number)[][],
  fieldType: ArrayFieldType
): { valid: boolean; errors: string[] } {
  const errors: string[] = [];
  const config = getFieldTypeConfig(fieldType);

  // Check all rows have correct dimensions
  arrayValue.forEach((row, idx) => {
    if (row.length !== config.dimensions) {
      errors.push(`Row ${idx + 1} has ${row.length} values, expected ${config.dimensions}`);
    }

    // For string types (like hex colors), validate format
    if (config.isStringType) {
      row.forEach((val, colIdx) => {
        if (typeof val !== 'string') {
          errors.push(`Row ${idx + 1}, ${config.columnLabels[colIdx]}: expected string, got ${typeof val}`);
        } else if (fieldType === 'color' && !isValidHexColor(val)) {
          errors.push(`Row ${idx + 1}, ${config.columnLabels[colIdx]}: invalid hex color "${val}"`);
        }
      });
    } else {
      // Check value ranges for numeric types
      if (config.valueRange) {
        row.forEach((val, colIdx) => {
          const numVal = val as number;
          if (numVal < config.valueRange![0] || numVal > config.valueRange![1]) {
            errors.push(
              `Row ${idx + 1}, ${config.columnLabels[colIdx]}: value ${numVal} outside range [${config.valueRange![0]}, ${config.valueRange![1]}]`
            );
          }
        });
      }
    }
  });

  return {
    valid: errors.length === 0,
    errors,
  };
}

/**
 * Validate hex color format
 */
function isValidHexColor(color: string): boolean {
  return /^#[0-9A-Fa-f]{6}$/.test(color);
}

/**
 * Create a new row with default values
 */
export function createDefaultRow(
  fieldType: ArrayFieldType,
  existingRows?: (string | number)[][]
): (string | number)[] {
  const config = getFieldTypeConfig(fieldType);

  // If there are existing rows, copy the last one
  if (existingRows && existingRows.length > 0) {
    return [...existingRows[existingRows.length - 1]];
  }

  // Otherwise create with default values
  return Array(config.dimensions).fill(config.defaultValue);
}

/**
 * Parse clipboard data (supports CSV, TSV, JSON)
 */
export function parseClipboardData(text: string): number[][] | null {
  try {
    // Try JSON first
    const parsed = JSON.parse(text);
    if (Array.isArray(parsed)) {
      return parsed;
    }
  } catch {
    // Not JSON, try CSV/TSV
    const lines = text.trim().split('\n');
    const result: number[][] = [];

    for (const line of lines) {
      // Try tab-separated, then comma-separated
      const values = line.includes('\t')
        ? line.split('\t')
        : line.split(',');

      const numbers = values.map(v => parseFloat(v.trim())).filter(n => !isNaN(n));

      if (numbers.length > 0) {
        result.push(numbers);
      }
    }

    if (result.length > 0) {
      return result;
    }
  }

  return null;
}

/**
 * Get default array value for a field type
 * Used when creating new geometries or editing without server data
 */
export function getDefaultArrayValue(fieldType: ArrayFieldType): (string | number)[][] {
  const config = getFieldTypeConfig(fieldType);

  const defaults: Record<ArrayFieldType, (string | number)[]> = {
    position: [0, 0, 0],
    direction: [0, 1, 0],
    color: ['#FF0000'], // Red hex color
    radius: [1.0],
    scale: [1.0],
    rotation: [0, 0, 0],
    size_2d: [1.0, 1.0], // Default width and height for Plane
    size_3d: [1.0, 1.0, 1.0], // Default width, height, depth for Box
    generic: [0],
  };

  const baseValue = defaults[fieldType];

  // Return single instance (one row) for 2D arrays
  return [baseValue];
}

/**
 * Get human-readable label for array shape
 * @param value - The array value (can be numbers or strings)
 * @param fieldType - The field type
 * @returns Human-readable label like "2 positions" or "1 box"
 */
export function getArrayShapeLabel(
  value: (string | number)[] | (string | number)[][],
  fieldType: ArrayFieldType
): string {
  const config = getFieldTypeConfig(fieldType);

  // Convert to 2D array format for consistent handling
  let rows: (string | number)[][];
  if (Array.isArray(value[0])) {
    rows = value as (string | number)[][];
  } else {
    // Flat array - chunk by dimensions
    const flatArray = value as (string | number)[];
    rows = [];
    for (let i = 0; i < flatArray.length; i += config.dimensions) {
      rows.push(flatArray.slice(i, i + config.dimensions));
    }
  }

  const count = rows.length;

  // Field-specific labels with singular/plural
  const labels: Record<ArrayFieldType, { singular: string; plural: string }> = {
    position: { singular: 'position', plural: 'positions' },
    direction: { singular: 'direction', plural: 'directions' },
    color: { singular: 'color', plural: 'colors' },
    radius: { singular: 'value', plural: 'values' },
    scale: { singular: 'value', plural: 'values' },
    rotation: { singular: 'rotation', plural: 'rotations' },
    size_2d: { singular: 'plane', plural: 'planes' },
    size_3d: { singular: 'box', plural: 'boxes' },
    generic: { singular: 'row', plural: 'rows' },
  };

  const label = labels[fieldType];
  const word = count === 1 ? label.singular : label.plural;

  return `${count} ${word}`;
}

/**
 * Get preview string for array tooltip
 * @param value - The array value (can be numbers or strings like hex colors)
 * @param maxItems - Maximum number of items to show
 * @returns Preview string like "[0, 0, 0], [1, 1, 1], ..." or hex color values
 */
export function getArrayPreview(
  value: (string | number)[] | (string | number)[][],
  maxItems: number = 2
): string {
  if (Array.isArray(value[0])) {
    const rows = value as (string | number)[][];
    const preview = rows.slice(0, maxItems).map(row =>
      `[${row.map(v => typeof v === 'number' ? v.toFixed(2) : v).join(', ')}]`
    ).join(', ');
    return rows.length > maxItems ? `${preview}, ...` : preview;
  } else {
    const flatArray = value as (string | number)[];
    const preview = flatArray.slice(0, maxItems * 3).map(v => typeof v === 'number' ? v.toFixed(2) : v).join(', ');
    return flatArray.length > maxItems * 3 ? `${preview}, ...` : preview;
  }
}
