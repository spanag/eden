

from .run_sim import runEden
# version.py gets auto-generated at build time
from .version import __version__, __version_info__
# afaik the only way to hide a "private" thing in Python is to hide it in a lambda, like with javascript, it's not a problem

# This works only with 'import *'
__all__ = []
__all__.extend([ 'runEden' ])
