import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import pandas as pd
import math

# custom functions and constants
from constants import *
from functions import *
from package4 import *

#Уеб Бъклинг
hf = zmax_front-zmin_front
hr = zmax_rear-zmin_rear
shearAvg = np.empty_like(y_linspace)
for i,y in enumerate(y_linspace):
    shearAvg[i] = V_func(y)/(hf*t+hr*t)

kv = 1.5 #<-------- chat

def webBuckling():
    return shearAvg*kv

#Скин Бъклинг
kc = 21
bplate = 37
def skinBuckling():
    return (((math.pi**2)*kc*Emod)/(12*(1-poission**2)))*((t/bplate)**2)

K = 32
#Колъмн Бъклинг
def kolumnBuckling():
    return K*(math.pi**2)*Emod*I/()
#Глобал Бъклинг