import matplotlib.pyplot as plt
import lib.plotting as fancyplots
import numpy as np
import pdb
import time
from matplotlib.widgets import Button
    #THIS IS RUN DURING AN ACTIVE FIGURE.

l = 10.0
class Index(object):
    def leendert(self):
        self.b = 5.0
        print('leendert?')
    def __init__(self):
        print(l)
        self.leendert()

    def square(self):
        return(self.a*self.a)


callback = Index()
print(callback.b)
