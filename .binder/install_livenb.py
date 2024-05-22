%%writefile .binder/install_livenb.py
import os, subprocess
from subprocess import STDOUT,PIPE
# basedir = '.'
def InstallCommonStuff(basedir):
	# !rm -rf sample_data; ls -A1 .* | xargs rm -rf # Clear all colab default stuff instead of clone or temp, if ever needed
	for x in [
		f'''apt-get update && sed 's/#.*//' "{basedir}/.binder/apt.txt" | DEBIAN_FRONTEND=noninteractive xargs apt-get install -y -qq''',
		f'pip install -r {basedir}/.binder/requirements.txt',
		f'rm -rf chrome-linux',
		f'HOME=$(pwd) bash {basedir}/.binder/postBuild'
	]:
		# Colab doesn't necessarily connect the subprocess output to console as it should... It does work on Jupyter though
		print(x)
		subprocess.call(['bash', '-c', x])

if 'COLAB_BACKEND_VERSION' in os.environ:
	InstallCommonStuff(basedir='.')
	# and Colab specific
	from google.colab import output
	output.enable_custom_widget_manager()
	# not sphinxcontrib i guess, exclude it? LATER

if 'DEEPNOTE_PROJECT_ID' in os.environ:
	InstallCommonStuff(basedir='../')
	# that's all

#TODO autodetect more cases that need installing stuff...
