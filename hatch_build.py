"""Custom build hook to build the frontend before packaging."""

import os
import shutil
import subprocess

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class FrontendBuildHook(BuildHookInterface):
    """Build the frontend using bun/npm before packaging."""

    PLUGIN_NAME = "frontend"

    def initialize(self, version: str, build_data: dict) -> None:
        """Build the frontend before the Python package is built.

        Note: The 'version' parameter is the build type ('standard' or 'editable'),
        NOT the package version. Use self.metadata.version for the actual version.
        """
        # Get actual package version from metadata
        pkg_version = self.metadata.version

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
            f"Building frontend with {pkg_manager} (version={pkg_version})..."
        )

        # Install dependencies
        subprocess.run(install_cmd, cwd=app_dir, check=True)

        # Build the frontend with version injected
        env = os.environ.copy()
        env["VITE_APP_VERSION"] = pkg_version
        # Increase Node.js heap size for memory-constrained environments (e.g., RTD)
        env["NODE_OPTIONS"] = env.get("NODE_OPTIONS", "") + " --max-old-space-size=4096"
        subprocess.run(build_cmd, cwd=app_dir, check=True, env=env)

        self.app.display_info("Frontend build complete.")
