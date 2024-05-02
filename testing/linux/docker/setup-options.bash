#!/bin/bash
# Some configuration options are placed here to be accessible from Docker using scripts.
# basically to allow mounting folders and read/write data as the same user, and also adapt to issues with rootless mode.

# NB: when rootless, bind mounts are always mapped to fake root, thus the following docker run needs to be run as "root" instead of current user https://github.com/moby/moby/issues/41497
# TODO move to general sudo_docker usage throughout...
SUDO_DOCKER="$(dirname "${BASH_SOURCE[0]}")/sudo_docker.bash"
AUTO_DOCKER_USER_OR_ROOT=$(id -u):$(id -g)

if $SUDO_DOCKER context show >/dev/null 2>&1 && [ "rootless" == $($SUDO_DOCKER context show) ]; then
    AUTO_DOCKER_USER_OR_ROOT=0:0
fi

# but also allow overriding
_DOCKER_USER_OR_ROOT=${DOCKER_USER_OR_ROOT:-AUTO_DOCKER_USER_OR_ROOT}
