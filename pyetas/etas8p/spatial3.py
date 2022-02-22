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

import numpy as np
from numba import jit


#%% ***************************************************************************


@jit(nopython=True)
def fr(r, w):
    D = w[0]
    q = w[1]
    return (1 - (1 + r**2 / D)**(1 - q)) / (2 * np.pi)

@jit(nopython=True)
def dD_fr(r, w):
    D = w[0]
    q = w[1]
    return (1 - q) * (1 + r**2 / D)**(-q) / D * r**2 / D / (2 * np.pi)

@jit(nopython=True)
def dq_fr(r, w):
    D = w[0]
    q = w[1]
    return (1 + r**2 / D)**(1 - q) * np.log(1 + r**2 / D) / (2 * np.pi)



