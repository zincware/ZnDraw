FROM python:3.12
SHELL ["/bin/bash", "--login", "-c"]

WORKDIR /usr/src/app
# required for h5py
RUN apt update && apt install -y gcc pkg-config libhdf5-dev build-essential

RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.0/install.sh | bash

RUN nvm install 20.9.0
RUN npm install -g bun

COPY ./ ./

RUN cd app && bun install && bun vite build && cd ..
RUN pip install -e .

EXPOSE 5003
# # # # The code to run when container is started:
ENTRYPOINT ["zndraw", "--port", "5003", "--no-browser"]
