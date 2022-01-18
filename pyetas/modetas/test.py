# -*- coding: utf-8 -*-
# pyetas
# Copyright (C) 2021-2022 Salvatore Iacoletti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
"""

import math
import numpy as np
from scipy import integrate
from pyetas.etas8p.lambdaf6 import (pdf_time_trunc, pdf_time_trunc_p, 
                                    pdf_time_trunc_c, integral_pdf_time_trunc,
                                    integral_pdf_time_trunc_p,
                                    integral_pdf_time_trunc_c)
from pdf_time_trunc2 import pdf_time_trunc2
from scipy import LowLevelCallable
from pyetas.etas8p.spatial5 import (fr, dq_fr, dgamma_fr, dD_fr,
                                    pdf_fr, dq_pdf_fr, dgamma_pdf_fr, dD_pdf_fr)


# LowLevelCallable.from_cython("__pycache__.pdf_time_trunc2.cpython-36", "pdf_time_trunc2")

ta = 1*365
c = 0.005
p = 1.1

import time
ti = time.time()
for _ in range(0,100):
    integral_pdf_time_trunc(0.0, c, p, ta, 0.0, 35.0)
print(time.time()-ti)

ti = time.time()
for _ in range(0,100):
    integrate.quad(pdf_time_trunc2, 0.0, 35.0, args=(c,p,ta))
print(time.time()-ti)

# from distutils.core import setup
# from Cython.Build import cythonize

# setup(ext_modules = cythonize('pdf_time_trunc2.pyx'))

gamma, D, q, mag = 0.8, 0.005, 2., 6.
def pdf_fr2(r, gamma, D, q, mag):
    out = (q-1) / (D*math.exp(gamma*mag)*math.pi) * np.power(1+r**2/(D*math.exp(gamma*mag)), -q)
    return out
integrate.quad(pdf_fr2, 0.0, 35.0, args=(gamma, D, q, mag))[0]*2*math.pi

import matplotlib.pyplot as plt
plt.figure()
plt.plot()





