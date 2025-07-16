// src/renderers/CustomRangeSlider.js
import React from 'react';
import { withJsonFormsControlProps } from '@jsonforms/react';
import { FormLabel, Slider, Box, Typography } from '@mui/material'; // Example with MUI for better UI

const CustomRangeSlider = ({ data, handleChange, path, label, schema, required }) => {

  const step = schema.step || 0.1;
  const value = data ?? schema.default ?? schema.minimum;

  return (
    <Box sx={{ marginBottom: 2 }}>
      <FormLabel required={required}>{label}</FormLabel>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, paddingX: 1 }}>
        <Slider
          value={value}
          onChange={(_event, newValue) => handleChange(path, newValue)}
          min={schema.minimum}
          max={schema.maximum}
          step={step}
          aria-labelledby="input-slider"
        />
        <Typography sx={{ minWidth: 35 }}>{value}</Typography>
      </Box>
    </Box>
  );
};

export default withJsonFormsControlProps(CustomRangeSlider);
