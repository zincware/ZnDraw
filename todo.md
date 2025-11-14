- [x] selection e.g. ctrl + a with multiple windows triggers infinite loop
- [x] left / right arrow should stop playback
- [ ] floor color is strange
- [x] shadow blur max 2
- [ ] settings on "empty" room not showing
- [x] add double click on bookmarks to edit
- [x] make bookmarks in python like `vis.bookmarks[0] = 10` just like the others.
- [x] vis.geometries[camera].curve = "curve" and vis.geometries[camera].curve_position = 0..1 to allow simulating camera along curve
- [x] have one camera per room, but by default no client is attached to it. Use `vis.geometries` for it, so multiple cameras are possible.
- [x] on new room set `uploading` blob to be true.
- [ ] delete room
- [ ] drag / drop, ctrl + v and ctrl + c and download / upload buttons
- [x] curve, have "default" as color
- [ ] default value for modifiers / selection / settings / analyses
- [x] link to github
- [x] SiMGen buttons and tutorial
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
- [x] new messages info should clear when leaving the room
- [x] include presenter lock into room metadata and lock:* endpoints, like metadata-lock.
- [ ] make the overview per user, remove "hidden" rooms but allow to restrain rooms for specific users. Remove user from url, but allow ?user=... to login. Have a anonymous user that only sees "public" rooms.
- [x] check adding geometries via python appear without reload (also in the geometries sidebar)
- [x] drag/drop not allowed when --file-browser not active
- [ ] make get_jwt_auth_headers a pytest fixture
- [x] split up routes.py -> search for all `@main.` decorator, make a list and organize into multiple files?
- [x] loading sidebar menus takes for ever!
- [x] can not change settings!
- [x] fix `vis.geometries` to only show keys and not full data
- [ ] being at frame 1 and then removing the first one might not update the frame to shift to 0. Check!
- [ ] instead of `apply_schema_feature` use `Field(json_schema_extra={"x-custom-type":"smiles-pack-array"})`
- [ ] path tracer does not accept changes to materials ( might only be with presets)
- [ ] drag and drop: ask for append, insert at cursor position (using virtual canvas, or so), or create a new room or even replace?
- [ ] admin / register global users
- [ ] auto-reconnect extensions
- [ ] disallow selection / hover on fix atoms constraints
- [ ] uploading triggers canvas re-renders (activate e.g. path tracer to see the effect)
- [ ] fix hosted url path in e.g. `copy python code`
- [ ] remove batch upload from read_file and filesystem. The `extend` has included batching!
- [x] consider where to compute connectivity: server / client / on-the-fly / delayed?
- [x] fix the presenter lock!!
- [ ] test register_fs, search `filesystem:load`
- [ ] let the celery workers directly write to the storage instead of using `vis.extend`. Add a way to emit?! Access redis / znsocket?!
- [ ] allow to invalidate all frames or just specific keys!
- [ ] write tests for append / extend / insert / delitem / setitem where you patch ZnDraw.lock to be nullcontext and assert lock errors!
- [x] lock should be per connection and not per user!
- [ ] update lock meta for editing data and for changing e.g. step (playback lock?), what else?
- [ ] use celery background thread instead of redis ttl and also check that locks are removed on socket connection close (required? short ttl will take care? locks independet of socket?)
- [ ] refactor workerId to use sessionId!
- [x] remove check_room_locked and other code duplication beyond the lock decorator!
- [ ] remove `f"room:{self.room_id}:locked"` has been replaced by fine-grained locks!
- [ ] check how to assign the SessionID to the flask socket connection! And clean alternative identifiers!! We have 
```py
join_room(f"room:{current_room}")
join_room(f"user:{user_name}")
```
which can not work, because one user can have multiple sessions in different rooms!
- [ ] unite locks im frontend! metaDataLock / presenterLock ... should be combined into one lock system!
- [ ] can not jump to frame when playing / playing mode is not exited when clicking on progress bar!
- [ ] Failed to apply adaptive resolution: 1 validation error for Sphere
resolution on remote filesystem uploads.
- [ ] make the lock Clickable? lock description only in tooltip? Make it look different? currently, the changing description messes up headbar layout!
- [ ] test register a room extension, which already exists as global extension! (UI!!) Make sure room/global in included, not just the name!
- [ ] test register a room filesystem, which already exists as global filesystem! (UI!!) Make sure room/global in included, not just the name!
- [ ] move "filesystem:register" to REST
- [ ] remove / unify `workerId` and `sessionId`
- [ ] Store tokens in HttpOnly cookies / use package for managing JWT token!
- [ ] add worker concurrency feature, currently set to 1!
- [ ] user a logs out and user b logs in, data might not be deleted properly!