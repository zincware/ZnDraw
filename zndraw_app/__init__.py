import eventlet

eventlet.monkey_patch()

import os
from engineio.payload import Payload

Payload.max_decode_packets = os.environ.get('ENGINEIO_PAYLOAD_MAX_DECODE_PACKETS', 16)
