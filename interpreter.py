#!/usr/bin/env python3

import code
import sys

from screen import *

console = code.InteractiveConsole(locals())
sys.exit(console.interact(banner='Screen interpreter'))
