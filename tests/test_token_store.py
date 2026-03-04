"""Tests for TokenStore and TokenEntry in server_manager."""

from __future__ import annotations

import json
import stat
from datetime import UTC, datetime

import pytest

from zndraw.server_manager import TokenEntry, TokenStore


@pytest.fixture
def token_store(tmp_path):
    """Create a TokenStore using a temporary directory."""
    return TokenStore(directory=tmp_path)


@pytest.fixture
def sample_entry():
    """Create a sample TokenEntry."""
    return TokenEntry(
        access_token="eyJ.test.token",
        email="user@example.com",
        stored_at=datetime(2026, 3, 1, 12, 0, 0, tzinfo=UTC),
    )


class TestTokenEntry:
    def test_fields(self, sample_entry):
        assert sample_entry.access_token == "eyJ.test.token"
        assert sample_entry.email == "user@example.com"
        assert sample_entry.stored_at == datetime(2026, 3, 1, 12, 0, 0, tzinfo=UTC)

    def test_serialization_roundtrip(self, sample_entry):
        data = sample_entry.model_dump_json()
        restored = TokenEntry.model_validate_json(data)
        assert restored == sample_entry


class TestTokenStoreGet:
    def test_nonexistent_server(self, token_store):
        assert token_store.get("http://localhost:8000") is None

    def test_no_tokens_file(self, token_store):
        """Returns None when tokens.json doesn't exist."""
        assert token_store.get("http://localhost:8000") is None

    def test_malformed_json_file(self, token_store):
        """Returns None gracefully when file contains invalid JSON."""
        tokens_file = token_store.directory / "tokens.json"
        tokens_file.write_text("not valid json{{{")
        assert token_store.get("http://localhost:8000") is None

    def test_malformed_entry(self, token_store):
        """Returns None when the entry doesn't match TokenEntry schema."""
        tokens_file = token_store.directory / "tokens.json"
        tokens_file.write_text(
            json.dumps({"http://localhost:8000": {"wrong_field": "value"}})
        )
        assert token_store.get("http://localhost:8000") is None


class TestTokenStoreSet:
    def test_set_and_get(self, token_store, sample_entry):
        token_store.set("http://localhost:8000", sample_entry)
        result = token_store.get("http://localhost:8000")
        assert result is not None
        assert result.access_token == sample_entry.access_token
        assert result.email == sample_entry.email
        assert result.stored_at == sample_entry.stored_at

    def test_overwrite_existing(self, token_store, sample_entry):
        token_store.set("http://localhost:8000", sample_entry)
        new_entry = TokenEntry(
            access_token="new.token",
            email="new@example.com",
            stored_at=datetime(2026, 3, 2, tzinfo=UTC),
        )
        token_store.set("http://localhost:8000", new_entry)
        result = token_store.get("http://localhost:8000")
        assert result is not None
        assert result.access_token == "new.token"
        assert result.email == "new@example.com"

    def test_multiple_servers(self, token_store, sample_entry):
        entry_a = sample_entry
        entry_b = TokenEntry(
            access_token="other.token",
            email="other@example.com",
            stored_at=datetime(2026, 3, 2, tzinfo=UTC),
        )
        token_store.set("http://localhost:8000", entry_a)
        token_store.set("https://remote.example.com", entry_b)

        result_a = token_store.get("http://localhost:8000")
        result_b = token_store.get("https://remote.example.com")

        assert result_a is not None
        assert result_a.access_token == "eyJ.test.token"
        assert result_b is not None
        assert result_b.access_token == "other.token"

    def test_file_permissions(self, token_store, sample_entry):
        """tokens.json should be owner-only read/write (0600)."""
        token_store.set("http://localhost:8000", sample_entry)
        tokens_file = token_store.directory / "tokens.json"
        mode = stat.S_IMODE(tokens_file.stat().st_mode)
        assert mode == 0o600


class TestTokenStoreDelete:
    def test_delete_existing(self, token_store, sample_entry):
        token_store.set("http://localhost:8000", sample_entry)
        token_store.delete("http://localhost:8000")
        assert token_store.get("http://localhost:8000") is None

    def test_delete_nonexistent(self, token_store):
        """Deleting a key that doesn't exist should not raise."""
        token_store.delete("http://localhost:9999")

    def test_delete_preserves_other_entries(self, token_store, sample_entry):
        other = TokenEntry(
            access_token="keep.me",
            email="keep@example.com",
            stored_at=datetime(2026, 3, 1, tzinfo=UTC),
        )
        token_store.set("http://localhost:8000", sample_entry)
        token_store.set("http://localhost:9000", other)
        token_store.delete("http://localhost:8000")

        assert token_store.get("http://localhost:8000") is None
        result = token_store.get("http://localhost:9000")
        assert result is not None
        assert result.access_token == "keep.me"
