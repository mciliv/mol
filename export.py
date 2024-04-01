import sys
from pathlib import Path
from util.googl import File

sys.path.append(str(Path('util')))
fs = File.at(Path('data/').absolute()).match()
