# Connection Dialog Snippet Fix

## Problem

The `ConnectionDialog` frontend component shows a Python snippet with `user="xxx@guest.user"` but no `password=`. Since the auth validation now requires both `user` + `password` or neither, copying and running this snippet raises:

> Missing --password (required when --user is provided)

## Design

### Single file change

`frontend/src/components/ConnectionDialog.tsx`

### Python snippet (always, no branching)

Remove `user=` from the snippet entirely. Without explicit credentials, the Python client's `resolve_token()` fallback chain handles auth:

1. Stored token from `~/.zndraw/tokens.json` (from `zndraw-cli auth login`)
2. Fresh guest session via `POST /v1/auth/guest`

```python
from zndraw import ZnDraw

vis = ZnDraw(
  url="{window.location.origin}/",
  room="{roomId}",
)
```

### Login hint (always shown)

Below the snippet, show the CLI login command with the current user's email:

> To authenticate as **{userName}**, run:
> ```bash
> zndraw-cli auth login --url {window.location.origin}/
> ```

### Session code (unchanged, with auth note)

Update the label to note that sessions require authentication:

> Access this browser session (requires authentication):

### No backend changes

The existing `resolve_token()` fallback chain and `zndraw-cli auth login` flow handle everything.

## Testing

- Playwright: open ConnectionDialog, verify snippet contains no `user=` parameter
- Playwright: verify login hint shows current user email and correct URL
- Playwright: verify session code section mentions authentication requirement
