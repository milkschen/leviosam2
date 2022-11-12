# Input version name
VER=${1}

if [[ -z ${VER} ]]; then
    echo "Error: input VER is not set"
    echo "Usage: sh update_docker.sh <version_name>"
    exit
fi

exit

# build the environment image
cd env
docker build --tag leviosam2_env .

# build the production image
cd ../production
docker build --tag leviosam2:${VER} --tag leviosam2:latest .

# push using the ${VER} tag
# docker tag leviosam2 naechyun/leviosam2:${VER}
docker push naechyun/leviosam2:${VER}

# push using the tag "latest"
# docker tag leviosam2 naechyun/leviosam2:latest
docker push naechyun/leviosam2:latest
