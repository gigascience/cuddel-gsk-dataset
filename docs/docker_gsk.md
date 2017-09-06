# Using Docker with the GSK data set 

## To run Docker server on Mac Pro

```bash
# Start Kitematic
# Run this command in a console
bash -c "clear && DOCKER_HOST=tcp://192.168.99.100:2376 DOCKER_CERT_PATH=/Users/peterli/.docker/machine/machines/default DOCKER_TLS_VERIFY=1 /bin/bash"

# Search for centos images
docker search centos

# Download Centos image
docker pull centos

# Run container using centos image and get shell access
docker run -it centos

# Commit changes to image instance
docker commit 31bb8717ba4e centos-test

# List running and inactive containers
docker ps -a

# Stop container before its removal
docker stop <CONTAINER ID>

# Remove container
docker rm centos

# After creating a dockerfile, you need to build it
docker build --rm --no-cache -t <IMAGE NAME> .

# Launch
docker run --privileged --name metax -v /sys/fs/cgroup:/sys/fs/cgroup:ro -p 80:80 -d metax

# Log into a running docker container
docker exec -it <CONTAINER ID> bash
```