- [x] selection e.g. ctrl + a with multiple windows triggers infinite loop
- [x] left / right arrow should stop playback
- [ ] floor color is strange
- [x] shadow blur max 2
- [ ] settings on "empty" room not showing
- [x] add double click on bookmarks to edit
- [x] make bookmarks in python like `vis.bookmarks[0] = 10` just like the others.
- [x] vis.geometries[camera].curve = "curve" and vis.geometries[camera].curve_position = 0..1 to allow simulating camera along curve
- [x] have one camera per room, but by default no client is attached to it. Use `vis.geometries` for it, so multiple cameras are possible.
- [ ] on new room set `uploading` blob to be true.
- [ ] delete room
- [ ] drag / drop, ctrl + v and ctrl + c and download / upload buttons
- [x] curve, have "default" as color
- [ ] default value for modifiers / selection / settings / analyses
- [ ] link to github
- [ ] SiMGen buttons and tutorial
- [x] sync mode default should be on
- [x] make the snackbar pop ups, e.g. `room locked` appear higher to not interfere with the progress bar.
- [ ] close zndraw button via "settings" menu upper right and different icons, must be distinguishable, also via name
- [ ] move to rest with invalidate and update invalidate to contain the data that changed, if small 
    - [x] bookmarks
    - [x] selections
    - [ ] frame selection
    - [ ] len?
    - [ ] queue?
- [x] blobs to plotyl
- [x] frames / atoms selection in ploty ( need to be tested )
- [ ] add prefetching and show prefetched frames in progress bar
- [x] `vis.geometries.pop(name)` to delete geometry only working on page reload
- [ ] allow "freezing" current configuration as refernce or show multiple configurations on top of each other
- [ ] `x-custom-type: dynamic-enum` support transformations via custom editor with changed to data like
- [ ] support delta storage, e.g. mapping is `1+5+11` and have something like `vis.update(1, "calc.energy", ...)` or so?
```json
{
    "source": "calc.forces",
}
// or with a transformation
{
  "transform": "norm", // predefined server side (or frontend?)
  "source": "calc.forces",
  "scale": 0.1,
  "clip": [0.01, 5.0]
}
// or more complex transformations (needed?)
{
  "transform": "compose",
  "ops": [
    {"transform": "norm", "source": "calc.forces"},
    {"transform": "scale", "factor": 0.1},
    {"transform": "clip", "min": 0.01, "max": 5.0}
  ]
}
```
- [ ] with geometries sidebar open, playback is much slower!
- [x] fix property inspector!
- [x] if there are now bonds / connectivity, still allow playback!
- [ ] loading a room with frames 0:10:100 should not set the room to be already loaded (could think about a clever way to reuse this data though)
- [x] selection no longer in sync (only on page refresh) also true for menus
- [ ] new messages info should clear when leaving the room
- [ ] include presenter lock into room metadata and lock:* endpoints, like metadata-lock.
- [ ] make the overview per user, remove "hidden" rooms but allow to restrain rooms for specific users. Remove user from url, but allow ?user=... to login. Have a anonymous user that only sees "public" rooms.
- [ ] check adding geometries via python appear without reload (also in the geometries sidebar)
- [x] drag/drop not allowed when --file-browser not active
- [ ] make get_jwt_auth_headers a pytest fixture
- [ ] split up routes.py -> search for all `@main.` decorator, make a list and organize into multiple files?
- [ ] loading sidebar menus takes for ever!
- [ ] can not change settings!
- [ ] fix `vis.geometries` to only show keys and not full data
- [ ] being at frame 1 and then removing the first one might not update the frame to shift to 0. Check!
- [ ] zarr storage issue, use `structures.h5` and then `rdkit2ase.pack` (2025-10-30 11:07:56,361 - zndraw.app.frame_routes - ERROR - Failed to write to Zarr store: Shape mismatch for key 'arrays.residuenames': existing shape (68,), new shape (320,).
Traceback (most recent call last):
  File "/work/fzills/uv-cache/archive-v0/zfPFBvtlCbskoEh1EZHjp/lib/python3.11/site-packages/zndraw/app/frame_routes.py", line 584, in append_frame
    storage.append(decoded_data)
  File "/work/fzills/uv-cache/archive-v0/zfPFBvtlCbskoEh1EZHjp/lib/python3.11/site-packages/zndraw/storage.py", line 468, in append
    self.extend([value])
  File "/work/fzills/uv-cache/archive-v0/zfPFBvtlCbskoEh1EZHjp/lib/python3.11/site-packages/zndraw/storage.py", line 471, in extend
    extend_zarr(self.group, values)
  File "/work/fzills/uv-cache/archive-v0/zfPFBvtlCbskoEh1EZHjp/lib/python3.11/site-packages/zndraw/storage.py", line 365, in extend_zarr
    _extend_recursive(root, entry, current_index, total_entries)
  File "/work/fzills/uv-cache/archive-v0/zfPFBvtlCbskoEh1EZHjp/lib/python3.11/site-packages/zndraw/storage.py", line 303, in _extend_recursive
    raise ValueError(
ValueError: Shape mismatch for key 'arrays.residuenames': existing shape (68,), new shape (320,).)
