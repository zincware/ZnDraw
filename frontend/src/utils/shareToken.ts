// In-memory, per-tab share-token store keyed by roomId.
// Never persisted (fresh share per tab lifetime).
const store = new Map<string, string>();

export function rememberShareToken(roomId: string, token: string): void {
  store.set(roomId, token);
}

export function getShareToken(roomId: string): string | undefined {
  return store.get(roomId);
}

export function parseShareFromLocation(roomId: string): string | undefined {
  const params = new URLSearchParams(window.location.search);
  const tok = params.get("share");
  if (tok) {
    rememberShareToken(roomId, tok);
    return tok;
  }
  return undefined;
}
