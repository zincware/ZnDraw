import dataclasses
import socket
import subprocess
import time


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
        self.process = subprocess.Popen(
            ["zndraw", self.file, "--port", str(self.port)], stdout=subprocess.PIPE
        )
        time.sleep(2)
        print(self.process.stdout.readline().decode("utf-8").strip())

    def _repr_html_(self):
        from IPython.display import IFrame

        return IFrame(
            src=f"http://127.0.0.1:{self.port}", width=self.width, height=self.height
        )._repr_html_()

    def __del__(self):
        self.process.kill()
