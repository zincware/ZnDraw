import { materialRenderers } from "@jsonforms/material-renderers";
import CustomColorPicker, {
	customColorPickerTester,
} from "../components/jsonforms-renderers/CustomColorPicker";
import CustomDynamicEnumWithColorPicker, {
	customDynamicEnumWithColorPickerTester,
} from "../components/jsonforms-renderers/CustomDynamicEnumWithColorPicker";
import CustomRangeSlider, {
	customRangeSliderTester,
} from "../components/jsonforms-renderers/CustomRangeSlider";
import CustomSmilesEditor, {
	customSmilesEditorTester,
} from "../components/jsonforms-renderers/CustomSmilesEditor";
import CustomSmilesPackEditor, {
	customSmilesPackEditorTester,
} from "../components/jsonforms-renderers/CustomSmilesPackEditor";
import DynamicEnumRenderer, {
	dynamicEnumTester,
} from "../components/jsonforms-renderers/DynamicEnumRenderer";
import MaterialEditor, {
	materialEditorTester,
} from "../components/jsonforms-renderers/MaterialEditor";
import OwnershipToggleRenderer, {
	ownershipToggleTester,
} from "../components/jsonforms-renderers/OwnershipToggleRenderer";
import PositionAttachmentRenderer, {
	positionAttachmentTester,
} from "../components/jsonforms-renderers/PositionAttachmentRenderer";
import PropertyInspectorRenderer, {
	propertyInspectorTester,
} from "../components/jsonforms-renderers/PropertyInspectorRenderer";
import Vertices2DRenderer, {
	vertices2DRendererTester,
} from "../components/jsonforms-renderers/Vertices2DRenderer";

/**
 * Custom JSONForms renderers including material renderers and custom components.
 * Order matters: higher priority renderers should come first.
 */
export const customRenderers = [
	...materialRenderers,
	{ tester: ownershipToggleTester, renderer: OwnershipToggleRenderer }, // Priority 10 - Ownership claim/release toggle
	{ tester: vertices2DRendererTester, renderer: Vertices2DRenderer }, // Priority 10 - 2D vertices editor for Shape
	{ tester: materialEditorTester, renderer: MaterialEditor }, // Priority 10 - Material editor with preset/object support
	{
		tester: positionAttachmentTester,
		renderer: PositionAttachmentRenderer,
	}, // Priority 10 - Position with optional CurveAttachment
	{ tester: dynamicEnumTester, renderer: DynamicEnumRenderer }, // Priority 10 - New unified renderer with transform support
	{ tester: propertyInspectorTester, renderer: PropertyInspectorRenderer }, // Priority 10 - Property Inspector
	{
		tester: customDynamicEnumWithColorPickerTester,
		renderer: CustomDynamicEnumWithColorPicker,
	}, // Priority 10 - Legacy renderer (will be removed)
	{ tester: customSmilesPackEditorTester, renderer: CustomSmilesPackEditor }, // Priority 6 - SMILES pack editor
	{ tester: customSmilesEditorTester, renderer: CustomSmilesEditor }, // Priority 5 - SMILES editor
	{ tester: customColorPickerTester, renderer: CustomColorPicker }, // Priority 5
	{ tester: customRangeSliderTester, renderer: CustomRangeSlider },
];

/**
 * Recursively traverses a JSON schema and injects dynamic enum values.
 * Only handles dynamic-geometries; dynamic-atom-props is handled by DynamicEnumRenderer.
 * @param schema The original JSON schema.
 * @param geometries Optional geometries object for dynamic-geometries feature.
 * @returns A new schema object with enums injected.
 */
export const injectDynamicEnums = (
	schema: any,
	geometries?: Record<string, any>,
): any => {
	// Create a deep copy to avoid mutating the original object from the react-query cache.
	const newSchema = JSON.parse(JSON.stringify(schema));

	const traverse = (obj: any) => {
		if (obj && typeof obj === "object") {
			// Inject geometry keys for fields with "dynamic-geometries" feature
			if (
				Array.isArray(obj["x-features"]) &&
				obj["x-features"].includes("dynamic-geometries") &&
				geometries
			) {
				const geometryFilter = obj["x-geometry-filter"];
				obj.enum = Object.keys(geometries).filter((key) => {
					if (!geometryFilter) return true;
					return geometries[key]?.type === geometryFilter;
				});
			}

			Object.keys(obj).forEach((key) => traverse(obj[key]));
		}
	};

	traverse(newSchema);
	return newSchema;
};
