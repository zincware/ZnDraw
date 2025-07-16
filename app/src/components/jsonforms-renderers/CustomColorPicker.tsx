// src/renderers/CustomColorPicker.js
import React from 'react';
import { withJsonFormsControlProps } from '@jsonforms/react';
import { FormLabel, Box } from '@mui/material'; // Example with MUI

const CustomColorPicker = ({ data, handleChange, path, label, schema, required }) => {
  const value = data ?? schema.default ?? '#000000';

  return (
    <Box sx={{ marginBottom: 2 }}>
      <FormLabel required={required}>{label}</FormLabel>
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, marginTop: 1 }}>
        <input
          type="color"
          value={value}
          onChange={(e) => handleChange(path, e.target.value)}
          style={{ width: '40px', height: '40px', border: 'none', background: 'none', cursor: 'pointer' }}
        />
        <code>{value}</code>
      </Box>
    </Box>
  );
};

export default withJsonFormsControlProps(CustomColorPicker);