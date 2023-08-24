import collections.abc
import dataclasses
from distributed import Client, Variable
from zndraw.data import DataHandler

import urllib.request


@dataclasses.dataclass
class ZnDraw(DataHandler):
    """The ZnDraw Interface"""
    client : Client = dataclasses.field(default_factory=Client)
    url: str = None

    def __post_init__(self):
        # TODO: start the server and get the url
        self.data_handler = DataHandler(self.client)
        if self.url is None:
            self.url = Variable("url").get()

    def display(self, index):
        """Display the atoms at the given index"""
        urllib.request.urlopen(self.url + f"/display/{index}")

if __name__ == "__main__":
    client = Client("tcp://127.0.0.1:57037")
    zndraw = ZnDraw(client=client)

    atoms = zndraw[0]

    data = [atoms.copy() for _ in range(100)]
    [x.rattle(0.1, seed=seed) for seed, x in enumerate(data)]

    zndraw[350:450] = data

    zndraw.display(405)

    # print(zndraw[400:])

    # # print(len(zndraw))
    # # print(zndraw)
    # # print(zndraw[0])

    # zndraw[405] = zndraw[1] # this is a dictionary so it can not be set to a value

    # zndraw[400:405] = zndraw[0:5] # this is a dictionary so it can not be set to a value
    # # # set request to zndraw.url / update
    # # print(urllib.request.urlopen(zndraw.url + "/update").read())

    # import pandas as pd

    # # Example DataFrame
    # data = {'atoms': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}
    # df = pd.DataFrame(data)

    # # Define the slice and values to be assigned
    # slicer = slice(2, 4)
    # values = [1, 2, 3]

    # # Perform the assignment using the .loc accessor and .values attribute
    # df.loc[slicer, "atoms"] = values

    # print(df)

