from celery import shared_task
import znh5md
from tqdm import tqdm
from zndraw.utils import update_colors_and_radii, atoms_to_dict
import logging

log = logging.getLogger(__name__)


@shared_task
def read_file() -> None:
    from zndraw import Client

    client = Client(room="testroom", url="http://localhost:5000")
    client.connect()
    io = znh5md.IO("/Users/fzills/tools/zndraw-communication-testing/structures.h5")
    frames = io[:]
    for atoms in tqdm(frames, desc="Uploading frames"):
        update_colors_and_radii(atoms)
        # atoms = atoms.repeat((4, 4, 4))
        client.append(atoms_to_dict(atoms))


@shared_task(bind=True)
def run_extension_task(self, room: str, category: str, extension: str, data: dict, user_id: str) -> dict:
    """Run a server-side extension task.

    This Celery task runs extensions that are provided by the server (provider="celery")
    rather than by client workers. It mimics the behavior of _on_task_run in the client.

    Args:
        self: Celery task instance (bound)
        room: The room ID where the task should run
        category: The extension category (e.g., 'modifiers', 'selections', 'settings')
        extension: The extension name
        data: The extension input data (schema-validated)
        user_id: The user ID who triggered the task

    Returns:
        dict: Task result with status and optional error message
    """
    from zndraw.zndraw import ZnDraw
    from flask import current_app

    log.info(f"Running server-side extension task: {category}/{extension} in room {room}")

    try:
        # Import the extension classes
        from zndraw.extensions.modifiers import modifiers
        from zndraw.extensions.selections import selections
        from zndraw.settings import settings

        category_map = {
            "selections": selections,
            "modifiers": modifiers,
            "settings": settings,
        }

        if category not in category_map:
            raise ValueError(f"Unknown category: {category}")

        if extension not in category_map[category]:
            raise ValueError(f"Unknown extension '{extension}' in category '{category}'")

        # Get the extension class and instantiate it with the provided data
        ext_class = category_map[category][extension]
        instance = ext_class(**data)

        # Create a ZnDraw client connected to the specific room
        # Use the server URL from the Flask app config
        vis = ZnDraw(room=room, url="http://localhost:5000", user=user_id)

        # Run the extension
        instance.run(vis)

        log.info(f"Successfully completed extension task: {category}/{extension}")
        return {"status": "success", "extension": extension, "category": category}

    except Exception as e:
        log.error(f"Error running extension task {category}/{extension}: {e}", exc_info=True)
        return {
            "status": "error",
            "message": str(e),
            "extension": extension,
            "category": category,
        }
