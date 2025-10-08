# Auto connect to an existing server on the same machine.
To reduce resources, the CLI should automatically connect to an existing server instance if one is available on the same machine.

Strategy
1. Check if ~/.zndraw/server.pid exists. If it does, read the PID and PORT and ZnDraw version from the file.
2. Verify if the process with the read PID is running.
3. Attempt to connect to the server on the read PORT.
4. If the connection is successful, use the existing server.
5. If any of the above checks fail, start a new server instance and write its PID and PORT and ZnDraw version to ~/.zndraw/server.pid.

# CLI options

- `zndraw file.xyz` default behavior: connect to existing server or start a new one if none exists.
- `zndraw --force-new-server file.xyz` always start a new server instance, ignoring any existing server.
- `zndraw --connect https://zndraw.myserver.com` connect to a specified remote server, bypassing local server checks.
- `zndraw --status` check if a local server is running and display its status without starting or connecting to it.
- `zndraw --shutdown` stop the local server if it is running, and remove the server.pid file.
- `zndraw --port 8080` specify a custom port to start a new server instance if no existing server is found.
