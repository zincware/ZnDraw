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
import { FrameMetadata } from "../myapi/client";

/**
 * Custom JSONForms renderers including material renderers and custom components.
 * Order matters: higher priority renderers should come first.
 */
export const customRenderers = [
  ...materialRenderers,
  {
    tester: customDynamicEnumWithColorPickerTester,
    renderer: CustomDynamicEnumWithColorPicker,
  }, // Priority 10
  { tester: customColorPickerTester, renderer: CustomColorPicker }, // Priority 5
  { tester: customRangeSliderTester, renderer: CustomRangeSlider },
];

/**
 * Recursively traverses a JSON schema and injects dynamic enum values.
 * @param schema The original JSON schema.
 * @param metadata The metadata object containing keys to inject.
 * @returns A new schema object with enums injected.
 */
export const injectDynamicEnums = (
  schema: any,
  metadata: FrameMetadata | undefined
): any => {
  // Create a deep copy to avoid mutating the original object from the react-query cache.
  const newSchema = JSON.parse(JSON.stringify(schema));

  const traverse = (obj: any) => {
    if (obj && typeof obj === "object") {
      // Check if the current object has our custom dynamic enum property
      if (obj["x-dynamic-enum"] === "AVAILABLE_ATOMS_KEYS" && metadata?.keys) {
        // Inject the keys from the metadata as enum values
        obj.enum = metadata.keys;
        // IMPORTANT: Preserve x-dynamic-enum and x-color-picker markers
        // They are needed by the custom renderer tester
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
