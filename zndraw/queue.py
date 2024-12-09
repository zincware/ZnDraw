import logging
import typing as t

import socketio.exceptions
import znsocket.exceptions

from zndraw.base import Extension

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


log = logging.getLogger(__name__)
TASK_RUNNING = "ZNDRAW TASK IS RUNNING"


def check_queue(vis: "ZnDraw") -> None:
    """Main loop to check and process modifier tasks for both private and public queues."""
    while True:
        if not vis._modifiers:
            vis.socket.sleep(1)
            continue
        try:
            process_modifier_queue(vis)
            process_public_queue(vis)
            vis.socket.sleep(1)
        except (znsocket.exceptions.ZnSocketError, socketio.exceptions.SocketIOError):
            log.warning("Connection to ZnDraw server lost. Reconnecting...")
            vis.socket.disconnect()
            vis.socket.sleep(1)


def process_modifier_queue(vis: "ZnDraw") -> None:
    """Process private modifier tasks in the queue."""
    modifier_queue = znsocket.Dict(
        r=vis.r,
        socket=vis._refresh_client,
        key=f"queue:{vis.token}:modifier",
    )

    for key in modifier_queue:
        if key in vis._modifiers:
            try:
                task = modifier_queue.pop(key)
                cls = vis._modifiers[key]["cls"]
                run_kwargs = vis._modifiers[key]["run_kwargs"]
                run_queued_task(vis, cls, task, modifier_queue, run_kwargs)
            except IndexError:
                pass


def process_public_queue(vis: "ZnDraw") -> None:
    """Process public modifier tasks in the public queue."""
    from zndraw import ZnDraw

    if not any(mod["public"] for mod in vis._modifiers.values()):
        return

    public_queue = znsocket.Dict(
        r=vis.r,
        socket=vis._refresh_client,
        key="queue:default:modifier",
    )

    for room, room_queue in public_queue.items():
        for key in room_queue:
            if key in vis._modifiers and vis._modifiers[key]["public"]:
                new_vis = ZnDraw(url=vis.url, token=room, r=vis.r)
                try:
                    task = room_queue.pop(key)
                    # run_queued_task(new_vis, key, task, room_queue)
                    cls = vis._modifiers[key]["cls"]
                    run_kwargs = vis._modifiers[key]["run_kwargs"]
                    run_queued_task(new_vis, cls, task, room_queue, run_kwargs)
                except IndexError:
                    pass
                finally:
                    new_vis.socket.sleep(1)
                    new_vis.socket.disconnect()


def run_queued_task(
    vis: "ZnDraw",
    cls: t.Type[Extension],
    task: dict,
    queue: znsocket.Dict,
    run_kwargs: dict | None = None,
) -> None:
    """Run a specific task and handle exceptions."""
    if not run_kwargs:
        run_kwargs = {}
    try:
        queue[TASK_RUNNING] = True
        cls(**task).run(vis, **run_kwargs)
    except Exception as err:
        vis.log(f"Error running `{cls}`: `{err}`")
    finally:
        queue.pop(TASK_RUNNING)
