import dataclasses
from zndraw.client import Client
from zndraw.extensions import Extension, ExtensionType
import typing as t
import requests

class _ExtensionStore(t.TypedDict):
    public: bool
    run_kwargs: dict|None
    extension: t.Type[Extension]


@dataclasses.dataclass
class ZnDraw:
    url: str
    room: str
    user: str

    def __post_init__(self):
        self.client = Client(url=self.url, room=self.room, user=self.user)
        self.client.connect()
        self._extensions: dict[str, _ExtensionStore] = {}

    @property
    def step(self) -> int:
        return self.client.step

    @step.setter
    def step(self, value: int):
        self.client.step = value

    def __len__(self) -> int:
        return len(self.client)

    @property
    def settings(self):
        return self.client.settings
    
    def register_extension(self, extension: t.Type[Extension], public: bool = False, run_kwargs: dict|None = None):
        # A WARNING ABOUT RUN_KWARGS!
        # If multiple workers are registering the same extension, work load will be distributed among them.
        # We do check that the extension schema are the same, but run_kwargs are not validated and can 
        # lead to unreproducible behavior if they are different among workers.
        if not hasattr(extension, "category"):
            raise ValueError("Extension must have a 'category' attribute and inherit from Extension base class.")
        name = extension.__name__
        if name in self._extensions:
            raise ValueError(f"Extension '{name}' is already registered.")
        if extension.category not in (cat.value for cat in ExtensionType):
            raise ValueError(f"Extension category '{extension.category}' is not valid. Must be one of {[cat.value for cat in ExtensionType]}.")
        self._extensions[name] = {
            "public": public,
            "run_kwargs": run_kwargs,
            "extension": extension
        }
        print(f"Registered extension '{name}' of category '{extension.category}'.")

        schema = extension.model_json_schema()
        if not public:
            response = self.client.sio.call("register:extension", {
                "name": name,
                "category": extension.category,
                "schema": schema,
                "public": False
            })
            print(f"Extension '{name}' registered with room '{self.room}'.")
        else:
            response = self.client.sio.call("register:extension", {
                "name": name,
                "category": extension.category,
                "schema": schema,
                "public": True
            })
            print(f"Extension '{name}' registered as public.")
        
        if response.get("status") != "success":
            raise RuntimeError(f"Failed to register extension '{name}': {response}")
        