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
def pdf_time_trunc2(t, c, p, ta):
    if t > ta:
        return 0.
    if not abs(p-1.) < 1e-06:
        pdf = (1-p)/((c+ta)**(1-p) - c**(1-p)) * (c+t)**(-p)
    else:
        pdf = (math.log(c+ta)-math.log(c))**(-1) * (c+t)**(-1)
    return pdf
