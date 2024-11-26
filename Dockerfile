FROM python:3.12
SHELL ["/bin/bash", "--login", "-c"]

WORKDIR /usr/src/app

# required for h5py
RUN apt update && apt install -y gcc pkg-config libhdf5-dev build-essential
RUN curl -fsSL https://bun.sh/install | bash

COPY ./ ./
RUN cd app && bun install && bun vite build && cd ..
RUN pip install -e .

ENTRYPOINT ["zndraw", "--port", "5003", "--no-browser"]
