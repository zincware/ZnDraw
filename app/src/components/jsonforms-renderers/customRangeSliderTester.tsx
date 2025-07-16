// src/testers/customRangeSliderTester.js
import { rankWith, schemaMatches } from '@jsonforms/core';

export const customRangeSliderTester = rankWith(
  5, // A high rank to ensure it overrides default number renderers
  schemaMatches((schema) => {
    return schema?.format === 'range' && schema.type === 'number';
  })
);
