# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""

import math
def pdf_time_trunc2(t, c, p, ta):
    if t > ta:
        return 0.
    if not abs(p-1.) < 1e-06:
        pdf = (1-p)/((c+ta)**(1-p) - c**(1-p)) * (c+t)**(-p)
    else:
        pdf = (math.log(c+ta)-math.log(c))**(-1) * (c+t)**(-1)
    return pdf
