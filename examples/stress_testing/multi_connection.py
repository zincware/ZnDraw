import typer 
import subprocess

app = typer.Typer()

@app.command()
def main(file: str, n: int, browser: bool = False):
    cmd = ["zndraw", file]
    if not browser:
        cmd.append("--no-browser")
    # run cmd n times in parallel
    for _ in range(n):
        subprocess.Popen(cmd)

if __name__ == "__main__":
    app()

