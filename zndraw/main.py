import dataclasses
import socket
import threading
import webbrowser

from zndraw import app, globals


@dataclasses.dataclass
class ZnDraw:
    file: str
    port: int = None
    width: int = 700
    height: int = 600

    def __post_init__(self):
        if self.port is None:
            sock = socket.socket()
            sock.bind(("", 0))
            self.port = sock.getsockname()[1]
            sock.close()

        globals.config.file = self.file
        self.thread = threading.Thread(target=app.run, kwargs={"port": self.port})
        self.thread.start()
        webbrowser.open(f"http://127.0.0.1:{self.port}")

    def _repr_html_(self):
        from IPython.display import IFrame

        return IFrame(
            src=f"http://127.0.0.1:{self.port}", width=self.width, height=self.height
        )._repr_html_()

    def __del__(self):
        self.thread.join()
