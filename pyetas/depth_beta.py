# -*- coding: utf-8 -*-
"""
@author: Salvatore
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
