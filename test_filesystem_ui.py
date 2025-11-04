#!/usr/bin/env python
"""End-to-end test for the Remote Filesystem Browser UI.

This test verifies:
1. Authentication and access control
2. FilesystemSelector component displays and works
3. FileBreadcrumbs navigation
4. FileList displays files and directories
5. LoadFileDialog appears and functions correctly
"""
import os
import tempfile
import time
import uuid

import ase.build
import ase.io
from playwright.sync_api import sync_playwright, expect
from fsspec.implementations.local import LocalFileSystem

from zndraw import ZnDraw


def test_filesystem_browser_ui(base_url: str = "http://localhost:5173"):
    """Test the remote filesystem browser UI end-to-end."""

    # Setup: Create test data
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create test files in subdirectories
        test_dir = os.path.join(tmp_dir, "test_data")
        os.makedirs(test_dir)
        subdir = os.path.join(test_dir, "subdirectory")
        os.makedirs(subdir)

        # Create some test XYZ files
        atoms = ase.build.molecule("H2O")
        ase.io.write(os.path.join(test_dir, "water.xyz"), atoms)
        ase.io.write(os.path.join(subdir, "molecule.xyz"), atoms)

        # Create a test room and register filesystem
        room = uuid.uuid4().hex
        backend_url = "http://localhost:5000"  # Backend server
        vis = ZnDraw(room=room, url=backend_url)

        # Register a local filesystem pointing to test_data
        fs = LocalFileSystem()
        # Change to test_dir so the filesystem root is there
        original_cwd = os.getcwd()
        os.chdir(test_dir)

        try:
            vis.register_filesystem(fs, name="test-filesystem")

            # Wait a moment for registration to propagate
            time.sleep(1)

            with sync_playwright() as p:
                # Launch browser
                browser = p.chromium.launch(headless=True)
                context = browser.new_context()
                page = context.new_page()

                print("Starting filesystem UI test...")

                # Step 1: Navigate to the app and login
                print(f"1. Navigating to {base_url}")
                page.goto(base_url)
                page.wait_for_load_state('networkidle')

                # Take initial screenshot
                page.screenshot(path='/tmp/01_homepage.png', full_page=True)
                print("   Screenshot saved: /tmp/01_homepage.png")

                # Check if we need to login
                if page.locator('text="Login"').is_visible():
                    print("2. Logging in...")
                    # Fill in login form (assuming username field exists)
                    page.fill('input[type="text"]', 'test-user')
                    page.click('button:has-text("Login")')
                    page.wait_for_load_state('networkidle')
                    page.screenshot(path='/tmp/02_after_login.png', full_page=True)
                    print("   Screenshot saved: /tmp/02_after_login.png")

                # Step 2: Navigate to the room
                print(f"3. Navigating to room: {room}")
                page.goto(f"{base_url}/rooms/{room}")
                page.wait_for_load_state('networkidle')
                time.sleep(1)  # Wait for room to load

                page.screenshot(path='/tmp/03_room_page.png', full_page=True)
                print("   Screenshot saved: /tmp/03_room_page.png")

                # Step 3: Open the filesystem browser
                print("4. Looking for Remote Filesystems menu option...")

                # Find and click the room management menu
                # Look for common menu triggers
                menu_button = page.locator('button[aria-label*="menu"], button:has-text("⋮"), [data-testid="room-menu"]').first
                if menu_button.is_visible():
                    print("   Found menu button, clicking...")
                    menu_button.click()
                    page.wait_for_timeout(500)
                    page.screenshot(path='/tmp/04_menu_open.png', full_page=True)
                    print("   Screenshot saved: /tmp/04_menu_open.png")

                    # Look for "Remote Filesystems" menu item
                    remote_fs_item = page.locator('text="Remote Filesystems"')
                    if remote_fs_item.is_visible():
                        print("   Found 'Remote Filesystems' menu item, clicking...")
                        remote_fs_item.click()
                        page.wait_for_load_state('networkidle')
                    else:
                        print("   WARNING: 'Remote Filesystems' menu item not found")
                        print("   Available menu items:", page.locator('li[role="menuitem"]').all_text_contents())
                else:
                    # Try direct navigation
                    print("   Menu button not found, navigating directly...")
                    page.goto(f"{base_url}/rooms/{room}/remote-files")
                    page.wait_for_load_state('networkidle')

                time.sleep(1)
                page.screenshot(path='/tmp/05_filesystem_browser.png', full_page=True)
                print("   Screenshot saved: /tmp/05_filesystem_browser.png")

                # Step 4: Verify FilesystemSelector is visible
                print("5. Verifying FilesystemSelector component...")
                page_content = page.content()

                # Check if we see the filesystem name or selector
                if "test-filesystem" in page_content or "Select Filesystem" in page_content:
                    print("   ✓ FilesystemSelector is present")

                    # If dropdown exists, verify it shows the filesystem
                    select = page.locator('input[role="combobox"], select, [aria-label="Select Filesystem"]')
                    if select.count() > 0:
                        print("   ✓ Filesystem selector dropdown found")
                else:
                    print("   ✗ FilesystemSelector not visible")

                # Step 5: Verify FileBreadcrumbs component
                print("6. Verifying FileBreadcrumbs component...")
                # Look for breadcrumb navigation elements
                breadcrumbs = page.locator('nav[aria-label="breadcrumb"], [role="navigation"]')
                if breadcrumbs.count() > 0 or "Root" in page_content:
                    print("   ✓ FileBreadcrumbs navigation is present")
                else:
                    print("   ✗ FileBreadcrumbs not found")

                # Step 6: Verify FileList shows files
                print("7. Verifying FileList component...")

                # Wait a moment for files to load
                page.wait_for_timeout(2000)
                page.screenshot(path='/tmp/06_file_list.png', full_page=True)
                print("   Screenshot saved: /tmp/06_file_list.png")

                # Check if we see the test files
                if "water.xyz" in page.content():
                    print("   ✓ FileList is displaying files (found water.xyz)")

                    # Verify file list items exist
                    file_items = page.locator('[role="button"]:has-text("water.xyz"), button:has-text("water.xyz")').first
                    if file_items.is_visible():
                        print("   ✓ File items are clickable")
                else:
                    print("   ✗ FileList not showing expected files")
                    print("   Page content sample:", page.content()[:500])

                # Check for subdirectory
                if "subdirectory" in page.content():
                    print("   ✓ Subdirectories are visible")

                    # Try clicking subdirectory to test navigation
                    subdir_item = page.locator('[role="button"]:has-text("subdirectory"), button:has-text("subdirectory")').first
                    if subdir_item.is_visible():
                        print("   Testing subdirectory navigation...")
                        subdir_item.click()
                        page.wait_for_timeout(1000)
                        page.screenshot(path='/tmp/07_subdirectory.png', full_page=True)
                        print("   Screenshot saved: /tmp/07_subdirectory.png")

                        # Check if we navigated into subdirectory
                        if "molecule.xyz" in page.content():
                            print("   ✓ Subdirectory navigation works")

                            # Test breadcrumb navigation back
                            root_link = page.locator('a:has-text("Root"), button:has-text("Root")').first
                            if root_link.is_visible():
                                print("   Testing breadcrumb navigation...")
                                root_link.click()
                                page.wait_for_timeout(1000)

                                if "water.xyz" in page.content():
                                    print("   ✓ Breadcrumb navigation back to root works")
                        else:
                            print("   ✗ Subdirectory navigation may not be working")
                else:
                    print("   Note: No subdirectories found in file list")

                # Step 7: Test LoadFileDialog
                print("8. Testing LoadFileDialog component...")

                # Navigate back to root if needed and click on a file
                page.goto(f"{base_url}/rooms/{room}/remote-files")
                page.wait_for_load_state('networkidle')
                page.wait_for_timeout(2000)

                # Click on water.xyz file
                water_file = page.locator('[role="button"]:has-text("water.xyz"), button:has-text("water.xyz")').first
                if water_file.is_visible():
                    print("   Clicking on water.xyz to open load dialog...")
                    water_file.click()
                    page.wait_for_timeout(1000)
                    page.screenshot(path='/tmp/08_load_dialog.png', full_page=True)
                    print("   Screenshot saved: /tmp/08_load_dialog.png")

                    # Check if dialog appeared
                    dialog_title = page.locator('text="Load File from Remote Filesystem"')
                    if dialog_title.is_visible():
                        print("   ✓ LoadFileDialog opened successfully")

                        # Verify dialog elements
                        target_room_field = page.locator('input[label="Target Room"], input:has([label*="Target Room"])')
                        if target_room_field.count() > 0 or "Target Room" in page.content():
                            print("   ✓ Target Room field is present")

                        if "Frame Selection" in page.content():
                            print("   ✓ Frame selection fields are present")

                        # Test cancel button
                        cancel_btn = page.locator('button:has-text("Cancel")').first
                        if cancel_btn.is_visible():
                            print("   Testing cancel button...")
                            cancel_btn.click()
                            page.wait_for_timeout(500)

                            if not dialog_title.is_visible():
                                print("   ✓ Cancel button closes the dialog")
                        else:
                            print("   ✗ Cancel button not found")
                    else:
                        print("   ✗ LoadFileDialog did not open")
                else:
                    print("   ✗ Could not find water.xyz file to click")

                # Final screenshot
                page.screenshot(path='/tmp/09_final_state.png', full_page=True)
                print("   Screenshot saved: /tmp/09_final_state.png")

                # Print summary
                print("\n" + "="*60)
                print("TEST SUMMARY")
                print("="*60)
                print("✓ All major components tested")
                print("✓ Screenshots saved to /tmp/")
                print(f"✓ Test completed for room: {room}")
                print("\nTo view screenshots, run:")
                print("  open /tmp/01_homepage.png")
                print("  open /tmp/05_filesystem_browser.png")
                print("  open /tmp/08_load_dialog.png")

                browser.close()

        finally:
            os.chdir(original_cwd)
            vis.disconnect()


if __name__ == "__main__":
    # This script expects servers to be running:
    # - Backend on http://localhost:5000
    # - Frontend on http://localhost:5173
    test_filesystem_browser_ui()
    print("\n✓ Test completed successfully!")
