import dataclasses

@dataclasses.dataclass
class JoinData:
    token: str
    uuid: str
    auth_token: str

@dataclasses.dataclass
class ModifierRunData:
    params: dict
    url: dict
    sid: str = None
    target: str = None

    @property
    def name(self) -> str:
        return self.params["method"]["discriminator"]

@dataclasses.dataclass
class AnalysisRunData:
    params: dict
    target: str = None