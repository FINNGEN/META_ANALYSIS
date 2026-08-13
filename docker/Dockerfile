FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8.2

ADD scripts/*.py /usr/local/bin/
ADD scripts/*R /usr/local/bin/

RUN chmod a+x /usr/local/bin/*.R && chmod a+x /usr/local/bin/*.py

COPY docker/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

COPY docker/install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R

RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/*

WORKDIR /
