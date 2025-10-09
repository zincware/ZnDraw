"""Test script for auto-discovery feature."""

from zndraw import ZnDraw

# Test 1: Auto-discovery (should work if server is running)
print("Test 1: Auto-discovery with url=None")
try:
    vis = ZnDraw(url=None, room="test_auto_discovery")
    print(f"✓ Successfully connected to: {vis.url}")
    print(f"  Room: {vis.room}")
except RuntimeError as e:
    print(f"✗ Expected error (no server running): {e}")

# Test 2: Explicit URL (should always work if server is at that URL)
print("\nTest 2: Explicit URL")
try:
    vis2 = ZnDraw(url="http://localhost:5000", room="test_explicit")
    print(f"✓ Successfully connected to: {vis2.url}")
    print(f"  Room: {vis2.room}")
except Exception as e:
    print(f"✗ Error: {e}")
