import collections.abc
import dataclasses
import urllib.request

from distributed import Client, Variable

from zndraw.data import DataHandler


@dataclasses.dataclass
class ZnDraw(DataHandler):
    """The ZnDraw Interface"""

    client: Client = dataclasses.field(default_factory=Client)
    url: str = None

    def __post_init__(self):
        # TODO: start the server and get the url
        self.data_handler = DataHandler(self.client)
        if self.url is None:
            self.url = Variable("url").get()

    def display(self, index):
        """Display the atoms at the given index"""
        urllib.request.urlopen(self.url + f"/cache/reset/")
        urllib.request.urlopen(self.url + f"/display/{index}")
        for idx in range(0, len(self), 100):
            urllib.request.urlopen(self.url + f"/cache/{idx}")


if __name__ == "__main__":
    client = Client("tcp://127.0.0.1:61587")
    zndraw = ZnDraw(client=client)

    atoms = zndraw[0]

    data = [atoms.copy() for _ in range(100)]
    [x.rattle(0.1, seed=seed) for seed, x in enumerate(data)]

    zndraw[350:450] = data

    zndraw.display(350)
