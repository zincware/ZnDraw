 TOP 3 CRITICAL ISSUES:

  1. Zustand Store Mass Re-renders (store.tsx)
    - Problem: 60+ state properties in flat object causing ALL components to re-render on ANY state change
    - Impact: Hundreds of unnecessary re-renders per second during playback
    - Fix: Use useShallow selectors instead of destructuring entire store
  // BAD - re-renders on every state change
  const { roomId, userName, selections } = useAppStore();

  // GOOD - only re-renders when these specific values change
  const { roomId, userName } = useAppStore(useShallow(state => ({
    roomId: state.roomId, userName: state.userName
  })));
  2. FigureWindow.tsx - 1137 line monster
    - Uncontrolled store subscriptions triggering expensive Plotly updates
    - DataGrid columns recreated every render with new function references
    - Needs to be split into smaller components with proper memoization
  3. Particles.tsx - useEffect dependency explosion
    - Single effect with 14+ dependencies doing expensive Float32Array operations
    - Component subscribes to 20+ store properties
    - Re-computes transforms on every render without memoization

  🟠 HIGH SEVERITY ISSUES:

  - SecondaryPanel.tsx:134-141 - Schema injection runs on every geometry update
  - SelectionsPanel.tsx:234-414 - DataGrid columns with inline renderCell functions cause 200+ cell
  re-renders
  - useSocketManager.ts - 20+ dependencies causing socket reconnections and race conditions
  - Canvas.tsx:112-128 - Entire Three.js scene destroyed/recreated on camera toggle (500-1000ms freeze)

  🟡 MEDIUM ISSUES:

  - ChatWindow.tsx - Markdown components object created fresh every render
  - App.tsx - QueryClient, theme, router recreated on every render
  - Missing React.memo on SelectionLayer, SecondaryPanel, SelectionsPanel, ChatWindow
  - Console.logs in production (SelectionLayer.tsx)

  Performance Impact Estimate:

  | Current State                       | After Fixes        |
  |-------------------------------------|--------------------|
  | ~30fps during playback              | ~60fps             |
  | 100-500+ unnecessary re-renders/sec | <50 re-renders/sec |
  | Janky Plotly updates                | Smooth transitions |
  | 500-1000ms camera toggle freeze     | <100ms             |

  Recommended Action Plan:

  Phase 1 (Immediate - 2-3 hours):
  1. Add useShallow to all useAppStore calls
  2. Wrap SelectionLayer, SecondaryPanel, ChatWindow in React.memo
  3. Remove console.logs from production

  Phase 2 (High Priority - 1 day):
  4. Refactor FigureWindow into smaller components
  5. Fix Particles.tsx useEffect complexity
  6. Stabilize SelectionsPanel DataGrid columns

  Phase 3 (Medium Priority - 2-3 days):
  7. Add useCallback to event handlers
  8. Refactor socket manager dependencies
  9. Fix Canvas camera toggle without unmounting

  Would you like me to start implementing fixes for any of these issues? I'd recommend starting with Phase 1
   - the Zustand store optimizations alone would give you a 70-90% reduction in unnecessary re-renders.w