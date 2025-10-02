- CLI upload files via ASE / H5MD, slicing
```bash
zndraw file.xyz --start 10 --stop 20 --step 2
zndraw file.xyz --index ':' # slice(None, None, None)
zndraw file.xyz --index '0:10:2' # slice(0, 10, 2)
zndraw file.xyz --index '7' # only the 7 configuration
zndraw file.xyz --index '-2' # only the second last configuration
zndraw file.xyz --index '0,2,4,6' # only configurations 0, 2, 4, 6
```
- streaming from file instead of from zarr store
```bash
zndraw file.xyz --stream
```
- follow a file (like tail -f) for live updates
```bash
zndraw file.xyz --follow
```
# template rooms
- new endpoints GET /api/templates 
```json
[
  {
    "id": "graphene-sheet",
    "name": "Graphene Sheet",
    "description": "A 5x5 periodic graphene supercell.",
    "thumbnail_url": "/api/rooms/graphene-sheet/thumbnail.png" // Optional but highly recommended!
  },
  {
    "id": "water-box",
    "name": "Box of Water",
    "description": "216 water molecules in a cubic box.",
    "thumbnail_url": "/api/rooms/water-box/thumbnail.png"
  }
]
```
and GET /api/rooms 
```json
[
  {
    "id": "room1",
    "template_id": "graphene-sheet", // Optional, if created from a template
  },
  {
    "id": "room2",
    "template_id": "water-box",
  }
]
```

and POST /api/rooms to create a new room from a template, if template_id is null, a blank room is created. If not given , the default template is used.
```json
{
  "name": "my-new-room-name",
  "template_id": "graphene-sheet"
}
```
with admin endpoints for
```
POST /api/rooms/{room_id}/promote # the room will be locked for editing for ever!
POST /api/templates/{template_id}/archive # the template will no longer be available for new rooms, but existing rooms will not be affected.
POST /api/templates/set-default/{template_id} # set the default template for new rooms
```


- upload to running server, if the server is not running on the same machine.
```bash
zndraw file.xyz --url http://myserver:1234
```
- detect running server on the same machine.
We are using a message broker, redis or celery, to communicate with the server. If a server is running on the same machine, we can detect it and use it as default. E.g. ~/.zndraw/config file.
If we want to run a different server and not the one on the same machine, specify --port will override this.
```bash
zndraw file.xyz --port 1234 # if 1234 is a running zndraw server, we will use it. If not, we will start a new server on port 1234.
```
```bash
zndraw file.xyz # if zndraw finds a running server on the same machine, it will use it, no matter the port.
# if not, it will start a new server on default port.
```

- multiple file uploads to different rooms. The first file will always be used as default template. UUIDs will be generated for all rooms.
```bash
zndraw file1.xyz file2.xyz file3.xyz
```
- multiple file uploads to the same room
```bash
zndraw file1.xyz file2.xyz file3.xyz --room myroom
```
- append to existing room
```bash
zndraw file1.xyz --room myroom
```

- jupyer support
- vs code support
- connectivity calculation (e.g. rdkit2ase?)
- znsocket / without redis
- in-memory zarr storage
- support / flatten? info / arrays / calc data
