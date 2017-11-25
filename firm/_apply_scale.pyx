'''
This turns out to be a bottleneck, so we pre-compile it with Cython.
'''

import numpy as np
cimport numpy as cnp

def _apply_scale(str seq, dict scale_values, cnp.ndarray output_buffer):
    
    cdef str aa
    cdef int i = 0
    cdef int n = len(seq)
    
    while i < n:
        
        aa = seq[i]
        
        if aa == '_':
            output_buffer[i] = np.nan
        else:
            output_buffer[i] = scale_values.get(aa, 0)
        
        i += 1