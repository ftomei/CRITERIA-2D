#!/bin/bash
docker stop watering_simulator
docker rm watering_simulator
docker build -t watering_simulator .
docker run --name watering_simulator --volume $(pwd):/home --detach -t watering_simulator
docker exec watering_simulator bash ./scripts/wrapper_experiments.sh $1 $2