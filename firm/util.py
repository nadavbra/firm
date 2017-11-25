from __future__ import absolute_import, division, print_function

import sys
import os
from datetime import datetime


### Project Functions ###

def log(message):
    print('FIRM|PID-%s [%s]: %s' % (os.getpid(), datetime.now(), message))
    sys.stdout.flush()
