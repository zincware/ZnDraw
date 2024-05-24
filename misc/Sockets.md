# Socket API of ZnDraw

- `room:frames` is used to request the frames from the specific room. It will use `default` room if the room does not exist.


## Main

### connect

- WebClient has a `session[token]` from routing through `\`, therefore it will
  be send to join the room `app.config["ROOM_HOSTS"][session[token]]`

### disconnect

### join(token)

- PyClient has no `session[token]`, therefore it provides a token to be sent to
  a room.

### exit

## Scene
