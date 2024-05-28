FROM continuumio/miniconda3
WORKDIR /usr/src/app
COPY ./ ./
RUN conda install python=3.11 nodejs
RUN conda install conda-forge::packmol

# Make RUN commands use the new environment:
SHELL ["conda", "run", "--no-capture-output", "-n", "base", "/bin/bash", "-c"]
RUN cd zndraw/static && npm install && cd ..
RUN pip install -e .[rdkit]


EXPOSE 5003
# # The code to run when container is started:
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "base", "zndraw", "--port", "5003", "--no-browser"]
