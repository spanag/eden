
# Sadly, Docker cannot really be run without super-user privileges

USAGE="$0: docker_command arg1 arg2 ..."

if [ "$#" == "0" ]; then
	echo "$USAGE"
	exit 1
fi

mydir=$(dirname "$0")
source $mydir/check_docker.bash
echo $maybe_sudo docker "$@"
${maybe_sudo} docker "$@"
