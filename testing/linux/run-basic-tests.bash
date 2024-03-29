#!/bin/bash
set -e

source "$(dirname "${BASH_SOURCE[0]}")/default-repo-path.bash"

echo "Testing Python module invocation..."
python3 -c "import eden_simulator; eden_simulator.runEden('$REPO_DIR/examples/LEMS_NML2_Ex25_MultiComp.xml')"

# test the console script as well
echo "Testing console script invocation..."
eden nml "$REPO_DIR/examples/LEMS_NML2_Ex25_MultiComp.xml"
