- [ ] selection e.g. ctrl + a with multiple windows triggers infinite loop
- [x] left / right arrow should stop playback
- [ ] floor color is strange
- [x] shadow blur max 2
- [ ] settings on "empty" room not showing
- [x] add double click on bookmarks to edit
- [ ] make bookmarks in python like `vis.bookmarks[0] = 10` just like the others.
- [ ] vis.geometries[camera].curve = "curve" and vis.geometries[camera].curve_position = 0..1 to allow simulating camera along curve
- [ ] have one camera per room, but by default no client is attached to it. Use `vis.geometries` for it, so multiple cameras are possible.
- [ ] on new room set `uploading` blob to be true.
- [ ] delete room
- [ ] drag / drop, ctrl + v and ctrl + c and download / upload buttons
- [ ] curve, have "default" as color
- [ ] default value for modifiers / selection / settings / analyses
- [ ] link to github
- [ ] SiMGen buttons and tutorial
- [ ] sync mode default should be on
- [ ] close zndraw button via "settings" menu upper right and different icons, must be distinguishable, also via name
- [ ] move to rest with invalidate and update invalidate to contain the data that changed, if small 
    - [ ] bookmarks
    - [ ] selections
    - [ ] frame selection
    - [ ] len?
    - [ ] queue?
- [ ] blobs to plotyl
- [ ] frames / atoms selection in ploty ( need to be tested )
- [ ] add prefetching and show prefetched frames in progress bar
- [ ] allow "freezing" current configuration as refernce or show multiple configurations on top of each other
- [ ] `x-custom-type: dynamic-enum` support transformations via custom editor with changed to data like
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