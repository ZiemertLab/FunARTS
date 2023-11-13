FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update
RUN apt-get install -y --no-install-recommends git wget unzip bzip2 build-essential ca-certificates pip gnupg default-jre-headless python3.8

# Install FunARTS Webapp

## Add related files
ADD . /funartswebapp/
RUN mkdir /results/
RUN mkdir /uploads/
RUN mkdir /reference

RUN cat > /run/funartswebapp.log

WORKDIR /funartswebapp/
RUN rm config/funartsapp_default.conf
RUN mv config/funartsapp_docker.conf config/funartsapp_default.conf

## Install pip depencies
RUN pip install -r requirements.txt
RUN apt-get install -y redis

RUN redis-server --daemonize yes

EXPOSE 5000
#CMD ["uwsgi","--ini","config/uwsgi_docker.conf"]

CMD redis-server --daemonize yes && uwsgi --ini config/uwsgi_docker.conf