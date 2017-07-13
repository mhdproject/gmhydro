 FROM python:2.7
 MAINTAINER garethcmurphy@gmail.com
 ENV PYTHONUNBUFFERED 1
 RUN mkdir /code
 WORKDIR /code
 COPY requirements.txt /code/
 RUN pip install -r requirements.txt
 COPY . /code/