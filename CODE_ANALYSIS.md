# Code Analysis Report

**Date:** 2025-10-16
**Status:** Ongoing maintenance and improvements

## Completed Items

The following code simplifications and architectural improvements have been completed:
- ✅ Three.js pooled object usage verified and used correctly
- ✅ Shared color processing extracted to `expandSharedColor()` utility
- ✅ Socket event handlers refactored with `createInvalidateHandler()` factory
- ✅ Removed unused `CustomToolbar` component from ArrayEditorDialog
- ✅ Removed unused processing functions (process3DData, process2DData, process1DData, processSizeAttribute)
- ✅ Python validators consolidated in BaseGeometry class
- ✅ Fixed opacity property access in Box.tsx

---

## Remaining Architecture Improvements

## Additional Considerations

### Cache Update Strategies

**Location:** `app/src/hooks/useSocketManager.ts`

**Issue:** Three different strategies for updating client state:

1. **Direct store updates:** `setSelections()`, `setBookmarks()`, etc.
2. **Cache invalidation:** `queryClient.invalidateQueries()`
3. **Direct cache updates:** `queryClient.setQueryData()`

**Example from `onChatMessageNew` (lines 177-194):**
```typescript
queryClient.setQueryData(["chat", roomId], (oldData: any) => {
  // Manually update cache
});
```

vs. `onSelectionsInvalidate` (lines 392-410):
```typescript
const response = await getAllSelections(roomId);
setSelections(response.selections); // Direct store update
```

**Recommendation:**
Always go:
User enteres data -> optimistically update store + invalidate cache + send to server -> fetch from server -> update store
- Prefer Cache invalidation with optimistic store updates for faster UX (not waiting for server to validate)
- Never use direct cache updates!
- Document the strategy for each data type

**Impact:** More predictable state management, easier debugging.

---

### 3.2 Geometry Invalidation Complexity

**Location:** `app/src/hooks/useSocketManager.ts:209-314`

**Issue:** 100+ line handler with nested conditionals, multiple fetch patterns, and complex state updates.

**Recommendation:** Split into focused handler functions (`handleGeometryDelete`, `handleGeometrySet`) to improve maintainability and testability.
