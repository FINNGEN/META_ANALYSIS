FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8.2

ADD scripts /META_ANALYSIS/scripts

RUN R -e "install.packages(c('openxlsx'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

WORKDIR /
