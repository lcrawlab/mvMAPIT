#!/usr/bin/env bash

CONTAINER_NAME=mvmapit

if [ ! "$(docker ps -q -f name=${CONTAINER_NAME})" ]; then
    if [ "$(docker ps -aq -f status=exited -f name=${CONTAINER_NAME})" ]; then
        docker rm ${CONTAINER_NAME}
    fi
    docker run --rm -ti \
    -e DISABLE_AUTH=true \
    -p 8787:8787 \
    --name ${CONTAINER_NAME} \
    mvmapit
fi