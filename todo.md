# Tier 1: Critical Performance & Stability Fixes
## Task 1: Create a Bulk setitem Endpoint on the Server.

Problem: The client's _simple_slice_assignment is dangerously inefficient. It makes a DELETE request for every old item and an INSERT request for every new item. vis[5:105] = new_data would result in 200 separate HTTP requests.

Solution: Create a new server-side endpoint, e.g., POST /api/rooms/<room_id>/frames/bulk-replace-slice. This endpoint should accept a start, stop, and a list of new frame data. It should perform the delete-and-insert logic atomically on the server in a single operation.

## Task 2: Optimize the Client's __setitem__ to Use the New Bulk Endpoint.

Problem: The client's _setitem_slice, _setitem_list, and _extended_slice_assignment methods are all inefficiently implemented with loops of single-item requests.

Solution:

Refactor _simple_slice_assignment to make a single call to the new bulk endpoint created in Task 1.

Refactor _setitem_list and _extended_slice_assignment to use a new bulk endpoint that accepts a list of indices and a corresponding list of data to replace. This can be a single POST request, which the server can then process efficiently.

## Task 3: Optimize the Server's Batch DELETE Logic.

Problem: The delete_frames_batch function loops through indices to be deleted and rebuilds the entire Redis sorted set on every single iteration. This is extremely inefficient for deleting large slices.

Solution: Refactor the function to perform the delete in one batch operation.

Get the current frame mapping from Redis once.

Calculate the final list of frames that should remain after the deletion.

Execute a single Redis transaction (pipeline): DELETE the old indices_key and ZADD the entire new mapping back in one go.

# Tier 2: High-Priority Enhancements
These tasks address major usability and deployment issues.

## Task 6: Improve Server-Side Logging.

Problem: The server code frequently uses print() for logging. This is unstructured and hard to manage in a production environment.

Solution: Replace all print() statements with the configured Flask logger (log.info, log.warning, log.error). This allows for proper log levels, formatting, and output redirection.

# Tier 3: Medium-Priority Refactoring & Optimization
These tasks will improve code quality and unlock further performance gains.

## Task 7: Optimize Bulk Frame Retrieval.

Problem: The get_frames endpoint iterates through the requested indices and loads each frame from Zarr individually.

Solution: Zarr supports "fancy indexing." Refactor the data retrieval to pass a list of physical indices to Zarr at once (storage.get([idx1, idx2, ...])), which can be significantly faster as it reads multiple chunks more efficiently.


## Task 12: Centralize Python Imports.

Problem: numpy is imported inside methods in the Client class.

Solution: Move all imports to the top of the file for clarity and to avoid repeated import overhead.

## Task 13: Add More Comprehensive Tests.

Uncomment and fix test_vis_slicing.py::test_delitem_list_with_duplicates and test_setitem_list_with_duplicates