FROM nvcr.io/nvidia/tensorflow:19.11-tf2-py3

LABEL maintainer="x25zeng@uwaterloo.ca"

# FOR PROD Copy PYTHON code into the container
# MUST run mvn clean install package first, otherwise the COPY ./target... stage will fail
COPY ./src/main/python /mstracer/src/main/python
COPY ./target/mstracer-1.0-SNAPSHOT-jar-with-dependencies.jar /mstracer/mstracer.jar
VOLUME /mstracer
# FOR PROD


# Compile Dependencies in order to run Pythonic Code
RUN apt-get update && \
	apt-get install -y maven && \
	pip install python-dotenv h5py tensorflow-gpu==2.1.0 scikit-learn==0.22.2.post1

# Export Python executables to the PATH
ENV PATH=/waterlooms/src/main/python:$PATH
ENV RUN_WITHIN_CONTAINER=TRUE

# Prevent asking for interactive status when installing the following package
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y install postgresql


# WORKDIR /waterlooms
#/usr/lib/postgresql/10/bin/pg_ctl -D /var/lib/postgresql/10/main >/waterlooms/db/logfile 2>&1 &

# Set the user for creation of waterlooms user
#USER postgres
#RUN /etc/init.d/postgresql start && \
#	psql --command "CREATE USER waterlooms WITH SUPERUSER PASSWORD 'waterlooms';" && \
#	createdb -O waterlooms waterlooms

# Revert to root for execution
USER root
