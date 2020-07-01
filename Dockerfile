#FROM gcr.io/finngen-refinery-dev/bioinformatics:0.5
FROM python:3.8

ADD . /META_ANALYSIS
RUN pip3 install -r /META_ANALYSIS/requirements/docker.txt
