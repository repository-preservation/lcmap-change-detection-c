FROM usgseros/ubuntu-python
MAINTAINER USGS LCMAP http://eros.usgs.gov

RUN apt-get install -y libgsl0-dev libgsl0ldbl gsl-bin \
    libmatio-dev libmatio2 gfortran

ENV CCDC_BIN /root/bin
RUN mkdir $CCDC_BIN

RUN git clone https://github.com/USGS-EROS/lcmap-change-detection-c.git
RUN cd lcmap-change-detection-c && \
    BIN=$CCDC_BIN make

ENTRYPOINT ["/root/bin/ccdc"]

