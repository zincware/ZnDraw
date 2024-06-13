FROM continuumio/miniconda3
WORKDIR /usr/src/app
COPY ./ ./
# required for h5py, chemfiles
RUN apt update && apt install -y gcc pkg-config libhdf5-dev build-essential
RUN conda install python=3.11 nodejs
# RUN conda install conda-forge::packmol
RUN npm install -g bun

# Make RUN commands use the new environment:
SHELL ["conda", "run", "--no-capture-output", "-n", "base", "/bin/bash", "-c"]
RUN cd app && bun install && bun vite build && cd ..
RUN pip install -e .[all]


EXPOSE 5003
# # The code to run when container is started:
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "base", "zndraw", "--port", "5003", "--no-browser"]
