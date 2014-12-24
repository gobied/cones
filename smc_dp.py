import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib as plt 
import pylab as pyl 
from math import *
import constants as consts
from matplotlib.backends.backend_pdf import PdfPages
from array import *

n_H = 0.3 #cm^-3
LAMBDA = 3550
SIGEXT = 9.599e-27 * n_H 
SIGSCA = 7.476e-27 * n_H 
SIGTOT = 7.518e-27 * n_H 
AVGCOS = 0.55799
ALBEDO = 0.77887

def readfile(fname):
    data = open(fname)
    theta = []
    diffcross = []
    for line in data:
        line = line.strip()
        line = line.split()
        line[0] = float(line[0])
        line[8] = float(line[8]) * n_H 
        theta.append(line[0])
        diffcross.append(line[8])
    theta = np.array(theta)
    diffcross = np.array(diffcross)
    return diffcross

def extrapolate(diffcross, th):
    while th > pi: 
        th = th - pi
    th = degrees(th)
    thint = int(th)
    thfrac = th - float(thint)
    return (thfrac * (diffcross[thint + 1] - diffcross[thint])) + diffcross[thint]
