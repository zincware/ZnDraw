"""RED PHASE TEST: Schema endpoint should return both public and private versions of the same extension.

This test is currently expected to FAIL because the backend only returns one version
per extension name due to using a dict with name as the key.

The fix requires the backend to use composite keys like "name:public" or return a list/dict
structure that can hold multiple variants of the same name.
"""

import pytest
import requests
from zndraw import ZnDraw
from zndraw.extensions import Extension, Category


class DuplicateNameExtension(Extension):
    """Test extension that will be registered in both public and private namespaces."""

    category = Category.MODIFIER
    test_value: int = 42

    def run(self, *args, **kwargs):
        """Dummy run method."""
        pass


def test_schema_returns_both_public_and_private_versions(server):
    """Schema endpoint should return BOTH public and private versions.

    This test registers the same extension in both namespaces and verifies that
    the schema endpoint returns both variants so the frontend can distinguish them.

    Expected behavior (after fix):
    - Schema should contain TWO entries for DuplicateNameExtension
    - One with public=true, one with public=false
    - Frontend can use composite keys to distinguish them

    Current behavior (before fix):
    - Schema contains only ONE entry for DuplicateNameExtension
    - Either public or private version (last registered wins)
    - Frontend cannot distinguish or select the other version
    """
    room = "test_dual_schema_room"

    # Connect and register BOTH versions
    vis = ZnDraw(
        url=server,
        room=room,
        user="test_user",
        auto_pickup_jobs=False
    )

    # Register private version first
    vis.register_extension(DuplicateNameExtension, public=False)

    # Register public version second
    vis.register_extension(DuplicateNameExtension, public=True)

    # Fetch schema from backend
    response = requests.get(
        f"{server}/api/rooms/{room}/schema/modifiers",
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"}
    )

    assert response.status_code == 200, f"Schema fetch failed: {response.text}"
    schemas = response.json()

    # Expected format: List of extension objects, just like filesystems
    # [
    #   {"name": "DuplicateNameExtension", "public": false, "schema": {...}, ...},
    #   {"name": "DuplicateNameExtension", "public": true, "schema": {...}, ...}
    # ]
    assert isinstance(schemas, list), (
        f"Expected schemas to be a list (like filesystems), but got {type(schemas).__name__}"
    )

    # Find all DuplicateNameExtension entries
    duplicate_extensions = [
        ext for ext in schemas
        if ext.get("name") == "DuplicateNameExtension"
    ]

    # Should have exactly 2 entries: one public, one private
    assert len(duplicate_extensions) == 2, (
        f"Expected 2 entries for DuplicateNameExtension (one public, one private), "
        f"but found {len(duplicate_extensions)}"
    )

    # Verify we have both public and private versions
    public_versions = [ext for ext in duplicate_extensions if ext["public"] is True]
    private_versions = [ext for ext in duplicate_extensions if ext["public"] is False]

    assert len(public_versions) == 1, (
        f"Expected exactly 1 public version, got {len(public_versions)}"
    )
    assert len(private_versions) == 1, (
        f"Expected exactly 1 private version, got {len(private_versions)}"
    )

    # Verify both have schemas
    assert "schema" in public_versions[0], "Public version should have schema"
    assert "schema" in private_versions[0], "Private version should have schema"

