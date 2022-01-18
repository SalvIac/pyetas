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
from scipy.stats import beta
import matplotlib.pyplot as plt


z = np.arange(0., 1.01, 0.01)

# Guo et al. 2015
eta = 10
m = 5
sei = 30
a = eta*m/sei+1
b = eta*(1-m/sei)+1
rv = beta(a, b)

fig, ax = plt.subplots(1, 1)
ax.plot(z, rv.pdf(z)/sei, 'k-', lw=2, label='frozen pdf')
# r = beta.rvs(a, b, size=1000)
plt.show()
