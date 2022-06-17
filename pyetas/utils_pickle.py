# -*- coding: utf-8 -*-
"""
@author: Salvatore Iacoletti
"""

import pickle

def load_pickle(name):
    with open(name + ".pkl", "rb") as f:
        return pickle.load(f)

def save_pickle(obj, name):
    with open(name + ".pkl", "wb") as f:
        pickle.dump(obj, file=f, protocol=pickle.HIGHEST_PROTOCOL)
