
# Sadly, Docker cannot really be run without sudo privileges

maybe_sudo=

# Check if docker can be run
docker --version
if [ $? -ne 0 ]; then
	echo "docker command not found; check PATH variable and/or installation"
	exit 1
fi

docker ps > /dev/null 2>&1
if [ $? -ne 0 ]; then
	echo "Docker will need sudo to run..."
	maybe_sudo=sudo
	# Run sudo to get access, once (not really required actually)
	# ${maybe_sudo} true 
fi
