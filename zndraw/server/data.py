import dataclasses

@dataclasses.dataclass
class JoinData:
    token: str
    uuid: str
    auth_token: str
