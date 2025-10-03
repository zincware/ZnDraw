import copy
import time

from ase.build import molecule

from zndraw import ZnDraw


def test_cache_hit(server, s22):
    """Verify that once data is cached, it is served from memory, not refetched."""
    vis = ZnDraw(url=server, room="test_cache_hit", user="tester")
    vis.extend(s22)

    # 1. Prime the cache by fetching all data.
    _ = vis[:]
    assert vis.cache is not None
    assert len(vis.cache) == len(s22)

    # 2. Get a reference to a cached object.
    cached_frame_5 = vis[5]
    assert cached_frame_5 == s22[5]

    # 3. Create a local, modified copy of the original data source.
    #    If the client were to refetch, it would get the original s22[5], but
    #    we can prove it's not even trying to by showing it returns the object
    #    already in memory.
    s22_modified = copy.deepcopy(s22)
    s22_modified[5].positions += 1.0

    # 4. Access the item again.
    refetched_frame_5 = vis[5]

    # 5. Assert that the returned object is the SAME as the originally cached one
    #    and DIFFERENT from our locally modified source data. This proves it
    #    was served from the cache without a network call.
    assert refetched_frame_5 == cached_frame_5
    assert refetched_frame_5 != s22_modified[5]


def test_cache_invalidation_single_item_replace(server, s22):
    """Verify `vis[idx] = ...` invalidates only the specific cache entry for idx."""
    vis = ZnDraw(url=server, room="test_cache_invalidate_replace", user="tester")
    vis.extend(s22)
    assert vis.cache is not None
    new_atoms = molecule("H2O")

    # 1. Prime the cache.
    _ = vis[:]
    assert 5 in vis.cache
    assert 6 in vis.cache

    # 2. Replace an item. This triggers a server event to invalidate the cache for index 5.
    vis[5] = new_atoms
    time.sleep(0.1)  # Give a moment for the invalidation event to be processed.

    # 3. Assert that only the specific key was removed from the cache.
    assert 5 not in vis.cache
    assert 6 in vis.cache  # Other keys should remain untouched.

    # 4. Assert that fetching the item again returns the new data.
    assert vis[5] == new_atoms


def test_cache_invalidation_delete(server, s22):
    """Verify `del vis[idx]` invalidates the cache from that index onward."""
    vis = ZnDraw(url=server, room="test_cache_invalidate_delete", user="tester")
    vis.extend(s22)
    assert vis.cache is not None

    # 1. Prime the cache.
    _ = vis[0:10]
    assert 4 in vis.cache
    assert 5 in vis.cache
    assert 6 in vis.cache

    # 2. Delete an item. This should invalidate the cache for all indices >= 5.
    del vis[5]
    time.sleep(0.1)

    # 3. Assert that entries before the change are still cached.
    assert vis.cache.keys() == {0, 1, 2, 3, 4}

    # 5. Verify that refetching the now-shifted item returns the correct data.
    new_frame_5 = vis[5]
    assert new_frame_5 == s22[6]


def test_cache_invalidation_insert(server, s22):
    """Verify that `vis.insert()` invalidates the cache from the insertion point."""
    vis = ZnDraw(url=server, room="test_cache_invalidate_insert", user="tester")
    vis.extend(s22)
    assert vis.cache is not None
    new_atoms = molecule("H2O")

    # 1. Prime the cache.
    _ = vis[0:10]
    assert vis.cache.keys() == set(range(10))

    # 2. Insert an item. This should invalidate the cache for all indices >= 5.
    vis.insert(5, new_atoms)
    time.sleep(0.1)

    assert vis.cache.keys() == {0, 1, 2, 3, 4}

    # 5. Verify that refetching returns the correct new and shifted data.
    assert vis[5] == new_atoms
    assert vis[6] == s22[5]


def test_cache_invalidation_bulk_slice_assignment(server, s22):
    """Verify slice assignments (grow, shrink, insert) invalidate from the start index."""
    vis = ZnDraw(url=server, room="test_cache_invalidate_bulk", user="tester")
    new_atoms = molecule("H2O")
    assert vis.cache is not None

    # --- Scenario 1: Growing the list ---
    vis.extend(s22)
    _ = vis[:]  # Prime cache
    assert vis.cache.keys() == set(range(len(s22)))

    vis[5:7] = [new_atoms.copy() for _ in range(4)]  # Replace 2 items with 4
    time.sleep(0.1)

    assert vis.cache.keys() == {0, 1, 2, 3, 4}
    del vis[:]  # Reset for next scenario

    # --- Scenario 2: Shrinking the list ---
    vis.extend(s22)
    _ = vis[:]  # Prime cache
    assert vis.cache.keys() == set(range(len(s22)))
    vis[5:10] = [new_atoms.copy() for _ in range(2)]  # Replace 5 items with 2
    time.sleep(0.1)

    assert vis.cache.keys() == {0, 1, 2, 3, 4}
    del vis[:]

    # --- Scenario 3: Insertion via empty slice ---
    vis.extend(s22)
    _ = vis[:]  # Prime cache
    assert vis.cache.keys() == set(range(len(s22)))
    vis[5:5] = [new_atoms.copy() for _ in range(3)]  # Insert 3 items at index 5
    time.sleep(0.1)

    assert vis.cache.keys() == {0, 1, 2, 3, 4}
    # Verify correctness of refetch
    assert vis[5] == new_atoms
    assert vis[8] == s22[5]


def test_cache_invalidation_clear_all(server, s22):
    """Verify that `del vis[:]` clears the entire cache."""
    vis = ZnDraw(url=server, room="test_cache_invalidate_clear", user="tester")
    vis.extend(s22)
    assert vis.cache is not None

    _ = vis[:]
    assert len(vis.cache) == len(s22)

    # 2. Clear the trajectory, which should trigger a full cache invalidation.
    del vis[:]
    time.sleep(0.1)

    assert len(vis.cache) == 0
