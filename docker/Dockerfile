FROM gcr.io/finngen-refinery-dev/bioinformatics:0.5

ADD META_ANALYSIS /META_ANALYSIS

WORKDIR /tmp

RUN wget https://github.com/atks/vt/archive/2187ff6347086e38f71bd9f8ca622cd7dcfbb40c.zip && \
    unzip 2187ff6347086e38f71bd9f8ca622cd7dcfbb40c.zip && \
    cd vt-2187ff6347086e38f71bd9f8ca622cd7dcfbb40c && make && \
    mv vt /usr/local/bin && cd .. && rm -rf vt-2187ff6347086e38f71bd9f8ca622cd7dcfbb40c

WORKDIR /
