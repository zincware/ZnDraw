import { Box, Skeleton } from "@mui/material";

/**
 * Skeleton loading placeholder for panels with schema selection dropdowns.
 * Used while initial schema data is being fetched.
 */
export const PanelSkeleton = () => (
  <>
    <Skeleton variant="rectangular" height={56} sx={{ mb: 2, borderRadius: 1 }} />
    <Skeleton variant="rectangular" height={40} sx={{ mb: 1 }} />
    <Skeleton variant="rectangular" height={40} sx={{ mb: 1 }} />
    <Skeleton variant="rectangular" height={40} sx={{ mb: 1 }} />
  </>
);

/**
 * Skeleton loading placeholder for JSONForms.
 * Used while form data or metadata is being fetched.
 */
export const FormSkeleton = () => (
  <Box sx={{ mt: 2 }}>
    <Skeleton variant="rectangular" height={56} sx={{ mb: 2, borderRadius: 1 }} />
    <Skeleton variant="rectangular" height={56} sx={{ mb: 2, borderRadius: 1 }} />
    <Skeleton variant="rectangular" height={56} sx={{ mb: 2, borderRadius: 1 }} />
    <Skeleton variant="text" width="60%" sx={{ mb: 1 }} />
    <Skeleton variant="rectangular" height={100} sx={{ borderRadius: 1 }} />
  </Box>
);

/**
 * Skeleton loading placeholder for selection tools accordion.
 * Specifically designed for the selections panel layout.
 */
export const SelectionToolsSkeleton = () => (
  <>
    <Skeleton variant="rectangular" height={48} sx={{ mb: 2, borderRadius: 1 }} />
    <Skeleton variant="rectangular" height={200} sx={{ mb: 2, borderRadius: 1 }} />
    <Skeleton variant="rectangular" height={48} sx={{ mb: 2, borderRadius: 1 }} />
  </>
);
