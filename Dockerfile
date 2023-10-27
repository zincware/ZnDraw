FROM continuumio/miniconda3
WORKDIR /usr/src/app
COPY ./ ./
RUN conda create -n myenv python=3.10 nodejs

# Make RUN commands use the new environment:
SHELL ["conda", "run", "--no-capture-output", "-n", "myenv", "/bin/bash", "-c"]
RUN cd zndraw/static && npm install && cd ..
RUN pip install .


EXPOSE 5003
# # The code to run when container is started:
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "myenv", "zndraw", "--port", "5003", "--no-browser", "--use-token"]
