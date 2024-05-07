FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8.2

ADD scripts/*.py /usr/local/bin/
ADD scripts/*R /usr/local/bin/

RUN chmod a+x /usr/local/bin/*.R && chmod a+x /usr/local/bin/*.py

RUN R -e "install.packages(c('openxlsx'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

WORKDIR /
