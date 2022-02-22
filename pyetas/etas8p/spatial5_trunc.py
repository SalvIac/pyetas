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


# PDF

def pdf_fr(r, gamma, D, q, mag):
    u = 0.95
    r_trunc = np.sqrt( D * np.exp(gamma*mag) * ((1-u)**(-1/(q-1))-1) )
    out = (q-1) / (D*np.exp(gamma*mag)*np.pi) * np.power(1+r**2/(D*np.exp(gamma*mag)), -q) / u
    out[r > r_trunc] = 0.
    return out


def dD_pdf_fr(r, gamma, D, q, mag):
    u = 0.95
    r_trunc = np.sqrt( D * np.exp(gamma*mag) * ((1-u)**(-1/(q-1))-1) )
    out = -((q - 1) * ((D + r**2 * np.exp(-gamma * mag))/D)**(-q - 1)*(D*np.exp(-gamma * mag) - \
            (q - 1) * r**2 * np.exp(-2 * gamma * mag)))/(np.pi * D**3) / u
    out[r > r_trunc] = 0.
    return out


def dq_pdf_fr(r, gamma, D, q, mag):
    u = 0.95
    r_trunc = np.sqrt( D * np.exp(gamma*mag) * ((1-u)**(-1/(q-1))-1) )
    sig = D*np.exp(gamma*mag)
    out = (np.power(1+r**2/sig, -q) * (1 - (q-1) * np.log(1+r**2/sig)))/(np.pi * sig) / u
    out[r > r_trunc] = 0.
    return out


def dgamma_pdf_fr(r, gamma, D, q, mag):
    u = 0.95
    r_trunc = np.sqrt( D * np.exp(gamma*mag) * ((1-u)**(-1/(q-1))-1) )
    out = -(mag * (q - 1) * ((D + r**2 * np.exp(-gamma * mag))/D)**(-q-1) * \
            (D * np.exp(-gamma * mag) - (q - 1) * r**2 * np.exp(-2 * gamma * mag)))/(np.pi * D**2) / u
    out[r > r_trunc] = 0.
    return out



#%%

if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    
    gamma, D, q, dm = 0.5, 0.001, 1.5, (7.5-3.5)
    
    # check 
    r = np.arange(0,2.,1e-4)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(r, pdf_fr(r, gamma, D, q, dm), label="fr")
    ax.plot(r, dgamma_pdf_fr(r, gamma, D, q, dm), label="dgamma_fr")
    ax.plot(r, dD_pdf_fr(r, gamma, D, q, dm), label="dD_fr")
    ax.plot(r, dq_pdf_fr(r, gamma, D, q, dm), label="dq_fr")
    ax.set_ylim([-80,80])
    plt.legend()
    plt.show()
    
    