from __future__ import absolute_import, division, print_function

import numpy as np

from .util import log

def apply_scale(seq, scale_values):

    if not isinstance(seq, str):
        seq = str(seq)
        
    output_buffer = np.empty(len(seq))
    _apply_scale(seq, scale_values, output_buffer)
    return output_buffer

def _pure_python_apply_scale(seq, scale_values, output_buffer):    
    for i, aa in enumerate(seq):
        if aa == '_':
            output_buffer[i] = np.nan
        else:
            output_buffer[i] = scale_values.get(aa, 0)

try:
    import pyximport; pyximport.install()
    from ._apply_scale import _apply_scale
except:
    log('apply_scale.py failed using Cython, will use a (slower) pure-python implementation instead.')
    _apply_scale = _pure_python_apply_scale
