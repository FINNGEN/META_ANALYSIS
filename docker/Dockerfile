FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.7

ADD scripts /META_ANALYSIS/scripts

RUN R -e "install.packages(c('rjson', 'ggpubr'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

WORKDIR /
