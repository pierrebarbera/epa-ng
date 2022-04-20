FROM ubuntu:20.04
RUN apt-get update && \
apt-get upgrade --assume-yes && \
DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC \
apt-get install --assume-yes autoconf automake autotools-dev bison build-essential cmake flex libtool zlib1g-dev && \
mkdir epa-ng
COPY CMakeLists.txt Makefile epa-ng/
COPY libs/ epa-ng/libs/
COPY src/ epa-ng/src/
COPY test/ epa-ng/test/
RUN cd epa-ng && make && cp bin/epa-ng /usr/local/bin/
CMD ["epa-ng", "--help"]
