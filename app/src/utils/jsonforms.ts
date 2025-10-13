import { materialRenderers } from "@jsonforms/material-renderers";
import CustomColorPicker, {
  customColorPickerTester,
} from "../components/jsonforms-renderers/CustomColorPicker";
import CustomRangeSlider, {
  customRangeSliderTester,
} from "../components/jsonforms-renderers/CustomRangeSlider";
import CustomDynamicEnumWithColorPicker, {
  customDynamicEnumWithColorPickerTester,
} from "../components/jsonforms-renderers/CustomDynamicEnumWithColorPicker";
import DynamicEnumRenderer, {
  dynamicEnumTester,
} from "../components/jsonforms-renderers/DynamicEnumRenderer";
import PropertyInspectorRenderer, {
  propertyInspectorTester,
} from "../components/jsonforms-renderers/PropertyInspectorRenderer";
import { FrameMetadata } from "../myapi/client";

/**
 * Custom JSONForms renderers including material renderers and custom components.
 * Order matters: higher priority renderers should come first.
 */
export const customRenderers = [
  ...materialRenderers,
  { tester: dynamicEnumTester, renderer: DynamicEnumRenderer }, // Priority 10 - New unified renderer
  { tester: propertyInspectorTester, renderer: PropertyInspectorRenderer }, // Priority 10 - Property Inspector
  {
    tester: customDynamicEnumWithColorPickerTester,
    renderer: CustomDynamicEnumWithColorPicker,
  }, // Priority 10 - Legacy renderer (will be removed)
  { tester: customColorPickerTester, renderer: CustomColorPicker }, // Priority 5
  { tester: customRangeSliderTester, renderer: CustomRangeSlider },
];

/**
 * Recursively traverses a JSON schema and injects dynamic enum values.
 * @param schema The original JSON schema.
 * @param metadata The metadata object containing keys to inject.
 * @param geometries Optional geometries object for dynamic-geometries feature.
 * @returns A new schema object with enums injected.
 */
export const injectDynamicEnums = (
  schema: any,
  metadata: FrameMetadata | undefined,
  geometries?: Record<string, any>
): any => {
  // Create a deep copy to avoid mutating the original object from the react-query cache.
  const newSchema = JSON.parse(JSON.stringify(schema));

  const traverse = (obj: any) => {
    if (obj && typeof obj === "object") {
      // NEW PATTERN: Check for x-custom-type="dynamic-enum" with "dynamic-atom-props" feature
      if (
        obj["x-custom-type"] === "dynamic-enum" &&
        Array.isArray(obj["x-features"]) &&
        obj["x-features"].includes("dynamic-atom-props") &&
        metadata?.keys
      ) {
        // Inject the keys from the metadata as enum values
        obj.enum = metadata.keys;
      }

      // NEW PATTERN: Check for x-custom-type="dynamic-enum" with "dynamic-geometries" feature
      if (
        obj["x-custom-type"] === "dynamic-enum" &&
        Array.isArray(obj["x-features"]) &&
        obj["x-features"].includes("dynamic-geometries") &&
        geometries
      ) {
        // Filter geometries by type if x-geometry-filter is present
        const geometryFilter = obj["x-geometry-filter"];
        const geometryKeys = Object.keys(geometries).filter((key) => {
          if (!geometryFilter) return true;
          const geometry = geometries[key];
          return geometry?.type === geometryFilter;
        });

        // Inject the filtered geometry keys as enum values
        obj.enum = geometryKeys;
      }

      // Continue traversing through the object's properties and array items
      Object.keys(obj).forEach((key) => {
        traverse(obj[key]);
      });
    }
  };

  traverse(newSchema);
  return newSchema;
};
