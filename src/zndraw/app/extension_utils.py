"""Utilities for working with extensions."""

import inspect

from zndraw.extensions import analysis, modifiers, selections


def get_server_extensions(category: str) -> set[str]:
    """Get list of server-side (Celery) extension names for a category.

    Server-side extensions are those defined in the zndraw.extensions modules.
    Client-side extensions are registered dynamically via the frontend.

    Args:
        category: Extension category ("modifiers", "selections", "analysis")

    Returns:
        Set of extension class names available on the server
    """
    module_map = {
        "modifiers": modifiers,
        "selections": selections,
        "analysis": analysis,
    }

    module = module_map.get(category)
    if not module:
        return set()

    # Get all classes defined in this module that have a 'run' method
    extension_classes = [
        cls
        for cls_name, cls in inspect.getmembers(module, inspect.isclass)
        if cls.__module__ == f"zndraw.extensions.{category}" and hasattr(cls, "run")
    ]

    return {cls.__name__ for cls in extension_classes}
