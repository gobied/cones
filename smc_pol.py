import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib as plt 
import pylab as pyl 
from math import *
import constants as consts
from matplotlib.backends.backend_pdf import PdfPages
from array import *


def readfile(fname):
    data = open(fname)
    theta = []
    pol = []
    for line in data:
        line = line.strip()
        line = line.split()
        line[0] = float(line[0])
        line[8] = float(line[8])
        theta.append(line[0])
        pol.append(line[8])
    theta = np.array(theta)
    pol = np.array(pol)
    return pol 
        
        
def extrapolate(pol, th):
    while th > pi:
        th = th - pi
    th = degrees(th)
    thint = int(th)
    thfrac = th - float(thint)
    return (thfrac * (pol[thint + 1] - pol[thint])) + pol[thint]
