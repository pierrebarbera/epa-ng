FROM ubuntu:18.04
RUN apt-get update
RUN apt-get -y install flex bison zlib1g-dev wget
RUN wget https://github.com/Pbdas/epa-ng/releases/download/v0.3.6/epa-ng .
RUN chmod u+x epa-ng
RUN mv epa-ng /usr/local/bin/

CMD epa-ng