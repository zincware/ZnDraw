# Deferred Items - Phase 01

## Pre-existing Test Failures

### test_close_clears_rooms (test_storage_asebytes.py)
- **Discovered during:** 01-02 Task 2
- **Status:** Pre-existing failure (also fails on prior commit without changes)
- **Issue:** `test_close_clears_rooms` expects `get_length()` to return 0 after `storage.close()`, but it returns 1
- **Impact:** Causes test pollution -- `test_router_extend_delegates_to_default` also fails when run after it due to a closed httpx client
- **Scope:** Out of scope for client package refactoring
