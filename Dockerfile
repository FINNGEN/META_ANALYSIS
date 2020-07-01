#FROM gcr.io/finngen-refinery-dev/bioinformatics:0.5
FROM alpine

pip3 install -r requirements/docker.txt
ADD META_ANALYSIS /META_ANALYSIS
