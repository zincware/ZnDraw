"""Custom build hook to build the frontend before packaging."""

import os
import shutil
import subprocess

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class FrontendBuildHook(BuildHookInterface):
    """Build the frontend using bun/npm before packaging."""

    PLUGIN_NAME = "frontend"

    def initialize(self, version: str, build_data: dict) -> None:
        """Build the frontend before the Python package is built."""
        # Get real version in editable mode
        if version == "editable":
            from hatch_vcs.version_source import VCSVersionSource

            version_source = VCSVersionSource(root=self.root, config={})
            version = version_source.get_version_data()["version"]

        app_dir = os.path.join(self.root, "app")

        if not os.path.isdir(app_dir):
            return

        # Determine which package manager to use
        if shutil.which("bun"):
            pkg_manager = "bun"
            install_cmd = [pkg_manager, "install"]
            build_cmd = [pkg_manager, "vite", "build"]
        elif shutil.which("npm"):
            pkg_manager = "npm"
            install_cmd = [pkg_manager, "install"]
            build_cmd = ["npx", "vite", "build"]
        else:
            raise RuntimeError("Neither bun nor npm found. Cannot build frontend.")

        self.app.display_info(
            f"Building frontend with {pkg_manager} (version={version})..."
        )

        # Install dependencies
        subprocess.run(install_cmd, cwd=app_dir, check=True)

        # Build the frontend with version injected
        env = os.environ.copy()
        env["VITE_APP_VERSION"] = version
        subprocess.run(build_cmd, cwd=app_dir, check=True, env=env)

        self.app.display_info("Frontend build complete.")
