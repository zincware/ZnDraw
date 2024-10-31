import typing as t


class Message(t.TypedDict):
    """A message to be sent to the client."""

    time: str
    msg: str
    origin: str | None
