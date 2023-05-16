#!/bin/bash
# Sadly, Docker cannot really be run without super-user privileges in the general case

USAGE="$0: docker_command arg1 arg2 ..."

if [ "$#" == "0" ]; then
	echo "$USAGE" >&2
	exit 1
fi

mydir=$(dirname "$0")
source $mydir/check_docker.bash
echo $maybe_sudo docker "$@" >&2
${maybe_sudo} docker "$@"
