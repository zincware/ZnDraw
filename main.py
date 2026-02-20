import uvicorn


def main() -> None:
    uvicorn.run(
        "zndraw.app:socket_app",
        host="0.0.0.0",
        port=8000,
        reload=True,
    )


if __name__ == "__main__":
    main()
