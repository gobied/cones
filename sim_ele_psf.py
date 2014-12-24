import numpy as np
import scipy
import scipy.signal
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.integrate import quad
import matplotlib as plt
import pylab as pyl
from math import *
import constants as consts
from matplotlib.backends.backend_pdf import PdfPages
from array import *
import glob
from optparse import OptionParser

parser = OptionParser()

parser.add_option('--plotInits', metavar = 'N', type = 'string', action = 'store', default = '', dest = 'plotInits', help='Program only plots the initial values')

parser.add_option('--fit', metavar = 'N', type = 'string',  action = 'store', default = '', dest = 'fit', help='Program fits model to data')

parser.add_option('--sim', action = 'store_true', default = False, dest = 'sim', help = 'Fit simultaneously the lateral and radial functions')

parser.add_option('--display', action = 'store_true', default = False, dest = 'display', help = 'Display the figures as they are plotted')

(options, args) = parser.parse_args()

drho = 0.001
dtheta = 0.001
dz = 0.001
EPSILON = 10**(-16)
VAR = 1.95697992754e-6
#d_L = 3828150
#norm = 1.26087e-15

def round_sigfigs(x, n):
    if x != 0:
        return round(x, int(n - ceil(log10(abs(x)))))
    else:
        return 0

def carttorho(x,y):
    return hypot(x, y)

def carttotheta(x,y):
    return atan2(y,x) - pi/2

def readfile(fname):
    data = open(fname)
    x = []
    y = []
    for line in data:
        line = line.strip()
        line = line.split()
        line[0] = float(line[0])
        line[1] = float(line[1])
        x.append(line[0])
        y.append(line[1])
    x = np.array(x)
    y = np.array(y)
    return x,y

def readinit(fname):
    inits = open(fname, 'r')
    ps= []
    lines = inits.readlines()
    for i in range(len(lines) - 1):
        ps.append(lines[i].strip().split())
    for i in range(len(ps)):
        for j in range(len(ps[1]) - 1):
            ps[i][j+1] = float(ps[i][j+1])

    return ps

def readsum(fname):
    summary = open(fname,'r')
    lines = summary.readlines()
    names = lines[0].strip().split()
    norms = lines[1].strip().split()[1:]
    #rlens = lines[2].strip().split()[1:]
    #plens = lines[3].strip().split()
    d_Ls = lines[5].strip().split()[1:]
    d_As = lines[6].strip().split()[1:]
    Ls = lines[8].strip().split()[1:]
    dists = lines[9].strip().split()[1:]
    for i in range(len(names)):
        norms[i] = float(norms[i]) * (4.7456221e3) * (1.8267533e-19)
        d_Ls[i] = float(d_Ls[i]) * 1000.
        d_As[i] = float(d_As[i])
        Ls[i] = float(Ls[i])
        dists[i] = float(dists[i])

    return names, norms, d_Ls, d_As, Ls, dists

def dA_perp(rho, theta, z):
    #Returns the surface area of a volume element projected onto the sphere
    A = cos(atan(z/rho)) * rho * dtheta * dz
    B = sin(atan(z/rho)) * rho * dtheta * drho
    return (A + abs(B))

def dV(rho, theta, z):
    #Returns the volume of a small volume element
    return rho * dtheta * drho * dz

def poltocart(rho, theta):
    #return a tuple x,y
    x0 = rho * cos(theta)
    y0 = rho * sin(theta)
    return x0, y0

def z_low(rho, theta, OMEGA, PHI):
    x0, y0 = poltocart(rho, theta)
    
    if OMEGA == PHI:
        a = 1 / (2 * x0 * tan(OMEGA))
        c = x0 / (tan(2 * OMEGA))
        return (a * y0 * y0 + c)

    if PHI == (pi - OMEGA):
        return (-1 * z_up(rho, theta, OMEGA, OMEGA))
   
    if PHI > OMEGA and PHI < (pi - OMEGA):
        anum = sqrt(2) * x0 * sin(OMEGA)
        aden = sqrt(cos(2 * OMEGA) - cos(2 * PHI))
        a = anum/aden
        OMEGA_app = atan(a / x0)
        #txtfile.write('theta = ' + str(degrees(theta)) + ' OMEGA = ' + str(degrees(OMEGA)) + ' OMEGA_app = ' + str(degrees(OMEGA_app)) + '\n')
        if theta > OMEGA_app or theta < -OMEGA_app:
            return 0
        bnum = x0 * sin(2 * OMEGA)
        bden = cos(2 * OMEGA) - cos(2 * PHI)
        b = bnum/bden
        zIII = x0 * sin(2 * PHI) /\
            (cos(2*OMEGA) - cos(2*PHI))
        return (zIII - ((b/a) * sqrt(a*a - y0*y0)))
    
    if PHI < OMEGA or (PHI > pi - OMEGA):
        if x0 >= 0:
            a = x0 / tan(PHI + OMEGA)
            b = x0 * cos(PHI + OMEGA) * sqrt(tan(OMEGA)) /\
                sqrt( abs( sin(2*PHI + OMEGA) * cos(OMEGA)))
            return (a * sqrt((y0*y0) + (b*b)) / b)
        elif x0 < 0:
            apr = abs(x0) / tan(PHI + OMEGA)
            bpr = abs(x0) * cos(PHI + OMEGA) * sqrt(tan(OMEGA)) /\
                sqrt(abs( sin(2*PHI + OMEGA) * cos(OMEGA)))
            c0 = ((abs(x0)/tan(OMEGA - PHI)) + (x0/tan(OMEGA+PHI)))/2
            return (c0 + (apr * sqrt((y0*y0) + (bpr*bpr)) / bpr))

def z_up(rho, theta, OMEGA, PHI):
    #Assumes cone extends to a maximum of sqrt(50) kpc
    x0, y0 = poltocart(rho, theta)
    
    if OMEGA == PHI:
        if (x0*x0 + y0*y0) > 100:
            return 0
        else:
            return sqrt(100 - x0 * x0 - y0 * y0)
    
    if PHI == (pi - OMEGA):
       return (-1 * z_low(rho, theta, OMEGA, OMEGA))

    if PHI > OMEGA and PHI < (pi - OMEGA):
        anum = sqrt(2) * x0 * sin(OMEGA)
        aden = sqrt(cos(2 * OMEGA) - cos(2 * PHI))
        a = anum/aden
        OMEGA_app = atan(a / x0)
        if theta > OMEGA_app or theta < -OMEGA_app:
            return 0
        bnum = x0 * sin(2 * OMEGA)
        bden = cos(2 * OMEGA) - cos(2 * PHI)
        b = bnum/bden
        zIII = x0 * sin(2 * PHI) /\
            (cos(2*OMEGA) - cos(2*PHI))
        return (zIII + ((b/a) * sqrt(a*a - y0*y0)))
    
    if PHI < OMEGA or (PHI > pi - OMEGA):
        if (100 - x0 * x0 - y0 * y0) < 0:
            return z_low(rho, theta, OMEGA, PHI)
        return sqrt((100 - x0 * x0 - y0 * y0))


def phase_thoms(rho, theta, z):
    #cons = (((consts.e * consts.e)/(consts.m_e * consts.c * consts.c))**2)/2
    cons = 2.82 * 10**(-13) 
    if rho == 0:
        th = 0
    elif z >= 0:
        th = (pi/2) - atan(z/rho)
    elif z < 0:
        th = (pi/2) + atan(z/rho)
    return cons*cons*(1 + cos(th)*cos(th))

def dflux_thoms(z, rho, theta, n0timesL, m, d_L, norm):
    '''
    The units are:
    A = erg/sec/(kpc)^3
    dV = kpc^3
    phase_thoms = cm^2
    hypot*hypot = kpc^2
    d_L = kpc
    df = erg * cm^2 / sec / kpc^4 = (1/(3.08567758*10^21))^4  erg/sec/cm^2
    '''

    A = n0timesL / ((hypot(z, rho)*hypot(z,rho))**m)
    phi = atan(z/rho)
    #B = abs(rho*dtheta*dz*cos(phi)) + abs(rho*dtheta*drho*sin(phi))
    #return (A * dV(rho, theta, z) * phase_thoms(rho, theta, z) /\
    #        (rho * dtheta * drho))
    df = (A * dV(rho, theta, z) * phase_thoms(rho, theta, z)) /\
        ( 4*pi * hypot(z, rho)*hypot(z, rho) * d_L*d_L )
    return df * ((1/(3.08567758*(10**21)))**4) / norm

def flux_thoms(rho, theta, n0timesL, OMEGA, PHI, m, d_L, norm):
    a = z_low(rho, theta, OMEGA, PHI)
    b = z_up(rho, theta, OMEGA, PHI)
    if a >= b: return 0
    result = 10**120 * quad(dflux_thoms, a, b, args=(rho, theta, n0timesL, m, d_L, norm))[0]
    return result

def rad_func(y, c, n0timesL, OMEGA, PHI, mthoms, d_L, norm, w):
    result = []
    #print (PHI > 23.) or (PHI < 20.) #or (mthoms > 1.5) or (mthoms < 0.) #or (n0timesL > 0.08) or (n0timesL < 0.04)
    if (OMEGA > pi/2) or (OMEGA < 0.) or (PHI < 0.) or (PHI > pi) or (mthoms > 5.) or (mthoms < 0.) or (PHI < OMEGA):

    #if (OMEGA > np.radians(16)) or (OMEGA < np.radians(13)) or (PHI > pi/2) or (PHI < -pi/2) or (mthoms > 1.5) or (mthoms < 0.) or (n0timesL > 0.08) or (n0timesL < 0.04):
    #if (OMEGA > pi/2) or (OMEGA < 0.) or (PHI < 0.) or (PHI > pi) or (mthoms > 5.) or (mthoms < 0.):

        for i in range(len(y)):
            result.append(1e10)
        return result

            
    yin = [-1.]
    d = y[1] - y[0]
    while yin[-1] < 1.:
        yin.append(yin[-1] + d)
    if len(yin)%2 == 0:
        yin.append(yin[-1] + d)

    xker, yker = np.meshgrid(np.linspace(-1,1,101), yin)
    
    ker = np.zeros(xker.shape)
    for i in range(xker.shape[0]):
        for j in range(xker.shape[1]):
            rr = sqrt( (xker[i,j] * xker[i,j]) + (yker[i,j] * yker[i,j]))
            if rr == 0.:
                ker[i,j] = 1.0
                continue
            ker[i,j] = _airy_func(rr, amplitude = 1.0, width = w)

    xx, yy = np.meshgrid(np.linspace(-1,1,101), y)
    zz = np.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            if carttorho(xx[i,j], yy[i,j]) == 0:
                print i, j
                continue
            zz[i,j] = flux_thoms(carttorho(xx[i,j], yy[i,j]),
                                 carttotheta(xx[i,j], yy[i,j]),
                                 n0timesL, OMEGA, PHI, mthoms, d_L, norm)


    ele = np.array([])
    for i in range(101):
        ele = np.append(ele, 0)
    pad = ele
    for i in range(((len(yin)-1)/2) - 1):
        pad = np.vstack((pad, ele))


    zz = np.vstack((pad, zz))
    zz = np.vstack((zz, pad))

    res = scipy.signal.fftconvolve(zz, ker, mode = 'valid')
    res = res.reshape(1,len(y))[0]
    res = np.array(res)

    return (res + c)

def perp_func(x, a, c, n0timesL, OMEGA, PHI, mthoms, d_L, norm, dist0, w):
    result = []
    #if (OMEGA > np.radians(60)) or (OMEGA < np.radians(50)) or (PHI > np.radians(40)) or (PHI < -pi/2) or (mthoms > 1.5) or (mthoms < 0.) or (n0timesL > 0.08) or (n0timesL < 0.04):
    if (OMEGA > pi/2) or (OMEGA < 0.) or (PHI < 0.) or (PHI > pi) or (mthoms > 5.) or (mthoms < 0.) or (PHI < OMEGA):

    #if (OMEGA > pi/2.) or (OMEGA < 0.) or (PHI < 0.) or (PHI > pi) or (a > 40.) or (a < 0.) or (mthoms > 5.) or (mthoms < 0.):
        for i in range(len(x)):
            result.append(1e10)
        return result
    

    x = x - a
    xx, yy = np.meshgrid(x, np.linspace(-1,1, 101))

    ker = np.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            rr = sqrt( (xx[i,j] * xx[i,j]) + (yy[i,j] * yy[i,j]))
            if rr == 0:
                ker[i,j] = 1.0
                continue
            ker[i,j] = _airy_func(rr, amplitude = 1.0, width = w)
    

    yy = yy + dist0
    zz = np.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            if carttorho(xx[i,j], yy[i,j]) == 0:
                print i, j
                continue
            zz[i,j] = flux_thoms(carttorho(xx[i,j], yy[i,j]),
                                 carttotheta(xx[i,j], yy[i,j]),
                                 n0timesL, OMEGA, PHI, mthoms, d_L, norm)
    
    ele = np.array([])
    if (len(x)%2) > 0:
        for i in range((len(x) - 1)/2):
            ele = np.append(ele, 0)
        pad = ele
        for i in range(100):
            pad = np.vstack((pad, ele))
        zz = np.hstack((pad, zz))
        zz = np.hstack((zz, pad))
    else:
        for i in range((len(x)/2) - 1):
            ele = np.append(ele, 0)
        eler = ele
        elel = np.append(ele, 0)
        padr = eler
        padl = elel
        for i in range(100):
            padl = np.vstack((padl, elel))
            padr = np.vstack((padr, eler))
        zz = np.hstack((padl, zz))
        zz = np.hstack((zz, padr))

    res = scipy.signal.fftconvolve(zz, ker, mode = 'valid')[0]

    return (res + c)

def make_perp(d_L, norm, dist0, w):
		def perp_func(x, a, c, n0timesL, OMEGA, PHI, mthoms):
			result = []
			#if (OMEGA > np.radians(16)) or (OMEGA < np.radians(13)) or (PHI > pi/2) or (PHI < -pi/2) or (mthoms > 1.5) or (mthoms < 0.) or (n0timesL > 0.08) or (n0timesL < 0.04):

			if (OMEGA > pi/2) or (OMEGA < 0.) or (PHI < 0.) or (PHI > pi) or (mthoms > 5.) or (mthoms < 0.) or (PHI < OMEGA):

			#if (OMEGA > pi/2.) or (OMEGA < 0.) or (PHI < 0.) or (PHI > pi) or (a > 40.) or (a < 0.) or (mthoms > 5.) or (mthoms < 0.):
				for i in range(len(x)):
					result.append(1e10)
				return result
			

			x = x - a
			xx, yy = np.meshgrid(x, np.linspace(-1,1, 101))

			ker = np.zeros(xx.shape)
			for i in range(xx.shape[0]):
				for j in range(xx.shape[1]):
					rr = sqrt( (xx[i,j] * xx[i,j]) + (yy[i,j] * yy[i,j]))
					if rr == 0:
						ker[i,j] = 1.0
						continue
					ker[i,j] = _airy_func(rr, amplitude = 1.0, width = w)
			

			yy = yy + dist0
			zz = np.zeros(xx.shape)
			for i in range(xx.shape[0]):
				for j in range(xx.shape[1]):
					if carttorho(xx[i,j], yy[i,j]) == 0:
						print i, j
						continue
					zz[i,j] = flux_thoms(carttorho(xx[i,j], yy[i,j]),
										 carttotheta(xx[i,j], yy[i,j]),
										 n0timesL, OMEGA, PHI, mthoms, d_L, norm)
			
			ele = np.array([])
			if (len(x)%2) > 0:
				for i in range((len(x) - 1)/2):
					ele = np.append(ele, 0)
				pad = ele
				for i in range(100):
					pad = np.vstack((pad, ele))
				zz = np.hstack((pad, zz))
				zz = np.hstack((zz, pad))
			else:
				for i in range((len(x)/2) - 1):
					ele = np.append(ele, 0)
				eler = ele
				elel = np.append(ele, 0)
				padr = eler
				padl = elel
				for i in range(100):
					padl = np.vstack((padl, elel))
					padr = np.vstack((padr, eler))
				zz = np.hstack((padl, zz))
				zz = np.hstack((zz, padr))

			res = scipy.signal.fftconvolve(zz, ker, mode = 'valid')[0]

			return (res + c)
		return perp_func

def make_rad(d_L, norm, w):
		def rad_func(y, c, n0timesL, OMEGA, PHI, mthoms):
			result = []
			#if (OMEGA > np.radians(16)) or (OMEGA < np.radians(13)) or (PHI > pi/2) or (PHI < -pi/2) or (mthoms > 1.5) or (mthoms < 0.) or (n0timesL > 0.08) or (n0timesL < 0.04):
			if (OMEGA > pi/2) or (OMEGA < 0.) or (PHI < 0.) or (PHI > pi) or (mthoms > 5.) or (mthoms < 0.) or (PHI < OMEGA):

			#if (OMEGA > pi/2) or (OMEGA < 0.) or (PHI < 0.) or (PHI > pi) or (mthoms > 5.) or (mthoms < 0.):
				for i in range(len(y)):
					result.append(1e10)
				return result

					
			yin = [-1.]
			d = y[1] - y[0]
			while yin[-1] < 1.:
				yin.append(yin[-1] + d)
			if len(yin)%2 == 0:
				yin.append(yin[-1] + d)

			xker, yker = np.meshgrid(np.linspace(-1,1,101), yin)
			
			ker = np.zeros(xker.shape)
			for i in range(xker.shape[0]):
				for j in range(xker.shape[1]):
					rr = sqrt( (xker[i,j] * xker[i,j]) + (yker[i,j] * yker[i,j]))
					if rr == 0.:
						ker[i,j] = 1.0
						continue
					ker[i,j] = _airy_func(rr, amplitude = 1.0, width = w)

			xx, yy = np.meshgrid(np.linspace(-1,1,101), y)
			zz = np.zeros(xx.shape)
			for i in range(xx.shape[0]):
				for j in range(xx.shape[1]):
					if carttorho(xx[i,j], yy[i,j]) == 0:
						print i, j
						continue
					zz[i,j] = flux_thoms(carttorho(xx[i,j], yy[i,j]),
										 carttotheta(xx[i,j], yy[i,j]),
										 n0timesL, OMEGA, PHI, mthoms, d_L, norm)


			ele = np.array([])
			for i in range(101):
				ele = np.append(ele, 0)
			pad = ele
			for i in range(((len(yin)-1)/2) - 1):
				pad = np.vstack((pad, ele))


			zz = np.vstack((pad, zz))
			zz = np.vstack((zz, pad))

			res = scipy.signal.fftconvolve(zz, ker, mode = 'valid')
			res = res.reshape(1,len(y))[0]
			res = np.array(res)

			return (res + c)

		return rad_func

def rad_gof(x, y, c, n0timesL, OMEGA, PHI, mthoms, d_L, norm, w):
    ex = rad_func(x, c, n0timesL, OMEGA, PHI, mthoms, d_L, norm, w)
    err = y - ex
    nu = len(y) - 5 - 1
    return sum((err*err)/(VAR*nu))

def perp_gof(x, y, a, c, n0timesL, OMEGA, PHI, mthoms, d_L, norm, dist, w):
    ex = perp_func(x, a, c, n0timesL, OMEGA, PHI, mthoms, d_L, norm, dist, w)
    err = y - ex
    nu = len(x) - 6 - 1
    return sum((err*err)/(VAR*nu))

def rad_err(x, y, c, n0timesL, OMEGA, PHI, mthoms, d_L, norm, w):
    ex = rad_func(x, c, n0timesL, OMEGA, PHI, mthoms, d_L, norm, w)
    err = ex - y
    return err

def perp_err(x, y, a, c, n0timesL, OMEGA, PHI, mthoms, d_L, norm, dist, w):
    ex = perp_func(x, a, c, n0timesL, OMEGA, PHI, mthoms, d_L, norm, dist, w)
    err = ex - y
    return err

def glob_err(p, x_r, x_p, y_r, y_p, d_L, norm, dist, w):
    c_r = p[0]
    a_p = p[1]
    c_p = p[2]
    n0timesL = p[3]
    OMEGA = p[4]
    PHI = p[5]
    mthoms = p[6]

    #OMEGA = np.radians(52)

    err_r = rad_err(x_r, y_r, c_r, n0timesL, OMEGA, PHI, mthoms, d_L, norm, w)
    err_p = perp_err(x_p, y_p, a_p, c_p, n0timesL, OMEGA, PHI, mthoms, d_L, norm, dist, w)

    return np.concatenate((err_r, err_p)) 

def _airy_func(rr, amplitude=1.0, width=.2):
    """
    For a simple radially symmetric airy function, returns the value at a given
    (normalized) radius
    """
    return amplitude * (2.0 * scipy.special.j1(2*1.61633*rr/width) / (2*1.61633*rr/width))**2

def main():
    radfiles = glob.glob('../dataCollection2/*rad*.dat')
    perpfiles = glob.glob('../dataCollection2/*perp*.dat')

    names, norms, d_Ls, d_As, Ls, dists = readsum('../dataCollection2/summary.txt')
    inits = readinit('inits.txt')

    tar = ['SDSSJ0149-0048', 'SDSSJ0210-1001', 'SDSSJ0319-0019', 'SDSSJ0319-0058', 'SDSSJ0321+0016', 'SDSSJ0759+1339', 'SDSSJ0841+2042', 'SDSSJ0842+3625', 'SDSSJ0858+4417', 'SDSSJ1039+4512', 'SDSSJ1040+4745']
    targs = {}
    for i in range(len(names)):
        targs[names[i]] = tar[i]


    if options.sim: prefix = 'sim'
    else: prefix = 'indep'
    suffix = ''
    if options.fit != '': suffix = '_fits'


    outf = open(prefix + suffix + '/fitparams.txt', 'w')

    properstring = True
    for i in range(len(options.fit)):
        properstring = properstring and (options.fit[i] in ',1234567890')
    #print properstring
    for i in range(len(options.plotInits)):
        properstring = properstring and (options.plotInits[i] in ',1234567890')
    #print properstring
    if options.fit == 'all' or options.plotInits == 'all':
        objs = range(11)
    elif properstring and options.fit != '':
        objs = np.array(options.fit.split(',')).astype(np.int) - 1
    elif properstring and options.plotInits != '':
        objs = np.array(options.plotInits.split(',')).astype(np.int) - 1
    else: 
        print 'Invalid option'
        return

    for i in objs:
        x_r, y_r = readfile(radfiles[i])
        x_r, y_r = np.delete(x_r, 0), np.delete(y_r, 0)

        x_p, y_p = readfile(perpfiles[i])

        name = inits[i][0]
        p = inits[i][1:]
        #p[4] = np.radians(p[4])
        #p[5] = np.radians(p[5])
        p = np.array(p)

        print '-----------------------------------------------------------'
        print name + ': ' + targs[name]
        print '-----------------------------------------------------------'
        print 'len(x_r) ', len(x_r)
        print 'len(x_p) ', len(x_p)

        c_r0, a_p0, c_p0, n0timesL0, OMEGA0, PHI0, mthoms0 = p

        w = 0.0367 * (1./3600.) * (pi / 180.) * 1000. * d_As[i] #kpc, 0.0367 arcsec is the min resolvable angular separation of Hubble

        if options.sim and options.fit != '':
            print 'Fitting both functions!'
            try:
                params = leastsq(glob_err, p, args = (x_r, x_p, y_r, y_p, d_Ls[i], norms[i], dists[i], w), maxfev = 1000)[0]
            except RuntimeError:
                print 'Possibly params not found. Rerun.'
                continue
            print 'Done fitting!'

            outf.write(name + '\t' + str(params) + '\n')# + name + '_p\t' + str(pparams) +'\n')

            #c_r, n0timesL_r, OMEGA_r, PHI_r, mthoms_r = rparams
            #a_p, c_p, n0timesL_p, OMEGA_p, PHI_p, mthoms_p = pparams
            c_r, a_p, c_p, n0timesL, OMEGA, PHI, mthoms = params
            print 'params= ', params
            #print 'pparams= ', pparams  
            n0timesL_r, n0timesL_p, OMEGA_r, OMEGA_p, PHI_r, PHI_p, mthoms_r, mthoms_p = n0timesL, n0timesL, OMEGA, OMEGA, PHI, PHI, mthoms, mthoms 

        #c_r, c_p, a_p, n0timesL, OMEGA, PHI, mthoms = params[0]

        if (not options.sim) and options.fit != '': 
            print 'Fitting the radial function!'
            try: rparams = curve_fit(make_rad(d_Ls[i], norms[i], w), x_r, y_r, [c_r0, n0timesL0, OMEGA0, PHI0, mthoms0], maxfev = 1000)[0]
            except RuntimeError:
                print 'Possibly rparams not found. Try again'
                continue
            print 'Fitting the lateral function!'
            try: pparams = curve_fit(make_perp(d_Ls[i], norms[i], dists[i], w), x_p, y_p, [a_p0, c_p0, n0timesL0, OMEGA0, PHI0, mthoms0], maxfev = 1000)[0]
            except RuntimeError:
                print 'Possibly pparams not found. Try again'
                continue
            print 'Done fitting!'

            outf.write(name + '_r\t' + str(rparams) + '\n' + name + '_p\t' + str(pparams) +'\n')

            c_r, n0timesL_r, OMEGA_r, PHI_r, mthoms_r = rparams
            a_p, c_p, n0timesL_p, OMEGA_p, PHI_p, mthoms_p = pparams
            print 'rparams= ', rparams
            print 'pparams= ', pparams  
    

        if options.plotInits != '':
            c_r, a_p, c_p, n0timesL, OMEGA, PHI, mthoms = c_r0, a_p0, c_p0, n0timesL0, OMEGA0, PHI0, mthoms0

            print 'c_r = ', c_r
            print 'a_p = ', a_p
            print 'c_p = ', c_p
            print 'n0timesL = ', n0timesL
            print 'OMEGA = ', OMEGA, degrees(OMEGA)
            print 'PHI = ', PHI, degrees(PHI)
            print 'mthoms = ', mthoms
    
            n0timesL_r, n0timesL_p, OMEGA_r, OMEGA_p, PHI_r, PHI_p, mthoms_r, mthoms_p = n0timesL, n0timesL, OMEGA, OMEGA, PHI, PHI, mthoms, mthoms 

        chi2_rad = rad_gof(x_r, y_r, c_r, n0timesL_r, OMEGA_r, PHI_r, mthoms_r, d_Ls[i], norms[i], w)
        chi2_perp = perp_gof(x_p, y_p, a_p, c_p, n0timesL_p, OMEGA_p, PHI_p, mthoms_p, d_Ls[i], norms[i], dists[i], w)
        print 'chi2_rad = ',chi2_rad
        print 'chi2_perp = ', chi2_perp
    
        x_rnew = np.linspace(x_r[0], x_r[-1], len(x_r))
        x_pnew = np.linspace(x_p[0], x_p[-1], len(x_p))
        y_rnew = rad_func(x_rnew, c_r, n0timesL_r, OMEGA_r, PHI_r, mthoms_r, d_Ls[i], norms[i], w)
        y_pnew = perp_func(x_pnew, a_p, c_p, n0timesL_p, OMEGA_p, PHI_p, mthoms_p, d_Ls[i], norms[i], dists[i], w)

        print 'len(x_rnew), len(y_rnew) = ', len(x_rnew), len(y_rnew)
        print 'len(x_pnew), len(y_pnew) = ', len(x_pnew), len(y_pnew) 

        #print y_rnew
        #print y_pnew

        plt.pyplot.subplot(121)
        plt.pyplot.plot( x_r, y_r, 'x', x_rnew, y_rnew, 'red')#, x_r, y_rnew2, 'green')
        plt.pyplot.ylim([-0.1, 1.0])
    
        plt.pyplot.xlabel('Distance/kpc')
        plt.pyplot.ylabel('Normlized Brightness')
   
        #plt.pyplot.text(24., 0.045, '$a = ' + str(round_sigfigs(a,3)) + '$')
        plt.pyplot.text(2., .95, name)#targs[name])
        plt.pyplot.text(x_rnew[-1] - x_rnew[-1]/5, .9, '$c_r = ' + str(round_sigfigs(c_r,3)) + '$')
        plt.pyplot.text(x_rnew[-1] - x_rnew[-1]/5, .85, '$n_0 L = ' + str(round_sigfigs(n0timesL_r,3)) + '$')
        plt.pyplot.text(x_rnew[-1] - x_rnew[-1]/5, 0.8, '$\Omega = '+str(round_sigfigs(degrees(OMEGA_r),3))+'^{\circ}$')
        plt.pyplot.text(x_rnew[-1] - x_rnew[-1]/5, 0.75, '$\Phi = ' + str(round_sigfigs(degrees(PHI_r), 3)) + '^{\circ}$')
        plt.pyplot.text(x_rnew[-1] - x_rnew[-1]/5, 0.7, '$m_{ele} = '+str(round_sigfigs(mthoms_r,3))+'$')
        plt.pyplot.text(x_rnew[-1] - x_rnew[-1]/5, 0.65, '$\chi ^2 = '+str(round_sigfigs(chi2_rad,3))+'$')
    
        plt.pyplot.subplot(122)
        plt.pyplot.plot(x_p, y_p, 'x', x_pnew, y_pnew, 'red')#, x_p, y_pnew2, 'green')
        plt.pyplot.ylim([min(y_p) -0.1*(max(y_p)-min(y_p)), max(y_p) + 0.05*max(y_p)])

        d = max(y_p) / 25.   
 
        plt.pyplot.xlabel('Distance/kpc')
        plt.pyplot.ylabel('Normalized Brightness')
 
        #plt.pyplot.text(2., .045, ) 
        plt.pyplot.text(x_pnew[-1] - x_pnew[-1]/5, max(y_p) + 0.05*max(y_p) - 3.*d, '$a_p = ' + str(round_sigfigs(a_p,3)) + '$')
        plt.pyplot.text(x_pnew[-1] - x_pnew[-1]/5, max(y_p) + 0.05*max(y_p) - 4.*d, '$c_p = ' + str(round_sigfigs(c_p,3)) + '$')
        plt.pyplot.text(x_pnew[-1] - x_pnew[-1]/5, max(y_p) + 0.05*max(y_p) - 5.*d, '$n_0 L = ' + str(round_sigfigs(n0timesL_p,3)) + '$')
        plt.pyplot.text(x_pnew[-1] - x_pnew[-1]/5, max(y_p) + 0.05*max(y_p) - 6.*d, '$\Omega = '+str(round_sigfigs(degrees(OMEGA_p),3))+'^{\circ}$')
        plt.pyplot.text(x_pnew[-1] - x_pnew[-1]/5, max(y_p) + 0.05*max(y_p) - 7.*d, '$\Phi = ' + str(round_sigfigs(degrees(PHI_p), 3)) + '^{\circ}$')
        plt.pyplot.text(x_pnew[-1] - x_pnew[-1]/5, max(y_p) + 0.05*max(y_p) - 8.*d, '$m_{ele} = '+str(round_sigfigs(mthoms_p,3))+'$')
        plt.pyplot.text(x_pnew[-1] - x_pnew[-1]/5, max(y_p) + 0.05*max(y_p) - 9.*d, '$\chi ^2 = '+str(round_sigfigs(chi2_perp,3))+'$')
    
        fig = plt.pyplot.gcf()
        fig.set_size_inches(16,6)
        plt.pyplot.savefig(prefix + suffix + '/' + prefix + '_' + name + suffix + '.png')
        if options.display:
            plt.pyplot.show(block = False)
            inp = raw_input("Hit Enter to Close")
        plt.pyplot.close()
    outf.close()

main()
