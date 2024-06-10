import eventlet.wsgi

eventlet.monkey_patch()  # MUST BE THERE FOR THE TESTS TO WORK

import pytest
import random
import socketio.exceptions

from zndraw.app import ZnDrawServer, FileIO

@pytest.fixture
def server():
    port = random.randint(10000, 20000)

    def start_server():
        fileio = FileIO()

        with ZnDrawServer(
        tutorial=None,
        auth_token=None,
        port=port,
        fileio=fileio,
        simgen=False,
        celery_worker=True,
        storage=None,
        ) as app:
            app.run(browser=False)

    thread = eventlet.spawn(start_server)

    # wait for the server to be ready
    for _ in range(100):
        try:
            with socketio.SimpleClient() as client:
                client.connect(f"http://localhost:{port}")
                client.disconnect()
                break
        except socketio.exceptions.ConnectionError:
            eventlet.sleep(0.1)
    else:
        raise TimeoutError("Server did not start in time")

    yield f"http://127.0.0.1:{port}"

    thread.kill()
