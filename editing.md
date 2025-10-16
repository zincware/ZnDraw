# Edit Mode
1. request meta lock
2. move spheres around -> others should see this!
3. release meta lock, this should first, save the new positions, etc. and then release the lock

For Particles, etc.. use https://tanstack.com/query/v5/docs/framework/react/guides/paginated-queries

Have a playback mode, that only goes to the next frame, once the current frame is fully loaded.

TODO: adaptive resolution. When reading, you can check how many particles, and reduce the resolution accordingly.
Progressbar in chat
Channel for specific frame / selection / ...

# Editing Mode for everything except "particles"
1. Only allow editing mode on objects where the positions are not dynamic, e.g. must be `number[][]`!
2. Enter editing mode and request lock (reject if not possible)
3. Create a virtual particle at the centroid of all selected object positions
4. Attach the transform controls to this virtual particle. OnChange, move update the matrix of the virtual particle, and update all selected objects accordingly. Have a debounce (like curve.tsx) and send the updated geometry to server (which will invalidate and let other clients fetch the new data)
5. leave editing mode and release lock

For particles, for now they need to be converted to `number[][]` first. We can have a button / modifer to create / append an atoms object from the geometry then. (consider having a Load into new geometry button to keep "particles").
