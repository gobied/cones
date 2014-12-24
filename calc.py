import math
import numpy as np
from scipy.integrate import quad
import glob

def dist_int(z, Om_m, Om_L):
	return ( 1/math.sqrt((Om_m * (1+z)**3) + Om_L) )

def round_sigfigs(x, n): 
    if x != 0:
        return round(x, int(n - math.ceil(math.log10(abs(x)))))
    else:
        return 0

def calcL(LOIII):
	M2500 = (LOIII + 0.62) / (-0.38)
	fnu = 10.**(-2. * (M2500 + 48.6) / 5.)
	Lnu = 4. * math.pi * 100. * 3.08567758e18 * 3.08567758e18 * fnu
	factor = (25./35.) ** (0.56)
	L = Lnu * factor * (3e10) / (2.5e-5)

	return L

def main():
	H0 = 100 * 0.678 * 1000 #(m/s)/Mpc
	Om_L = 0.692
	Om_m = 0.308

	z = []
	LOIII = []
	L = []
	mval = []
	rpix = []
	ppix = []
	dL = []
	dA = []
	dists = []

	summary = open('summary.txt', 'w')

	files = glob.glob('/Users/georgesobied/Documents/University/AstroResearch/fitssect/*.dat')
	print files

	for f in files:
		if f.find('perp') >= 0:
			inf = open(f, 'r')
			lines = inf.readlines()
			pix = float(lines[-1].strip().split()[0])

			ppix.append(pix)


		elif f.find('rad') >= 0:
			inf = open(f, 'r')
			lines = inf.readlines()
			val = float(lines[3].strip().split()[1])
			pix = float(lines[-1].strip().split()[0])

			mval.append(val)
			rpix.append(pix)


	print 'mval ', mval
	print 'ppix ', ppix
	print 'rpix ', rpix
	
	summary.write('\t\t')
	for j in range(len(mval)):
		name = files[2 * j][ files[2*j].rfind('/') + 1: -9 ].ljust(11)
		summary.write(name + '\t')
	summary.write('\n mval\t')
	for j in range(len(mval)):
		summary.write(str(mval[j]).ljust(11) + '\t')
	summary.write('\n ppix\t')
	for j in range(len(ppix)):
		summary.write(str(ppix[j]).ljust(11) + '\t')
	summary.write('\n rpix\t')
	for j in range(len(rpix)):
		summary.write(str(rpix[j]).ljust(11) + '\t')

	mval = np.array(mval)
	ppix = np.array(ppix)
	rpix = np.array(rpix)

	i = 0
	for f in files:
		name = f[f.rfind('/') + 1:]
		rawf = open(f,'r')
		outf = open('dataCollection/' + name[:-4] + '_norm.dat', 'w')
	
		if f.find('perp') >= 0:
			ztemp = float(rawf.readline().strip().split()[1])
			LOIIItemp = float(rawf.readline().strip().split()[1])
			lentemp = float(rawf.readline().strip().split()[1])
			disttemp = float(rawf.readline().strip().split()[1])
			dLtemp = 3e8 * (1 + ztemp) * quad(dist_int, 0, ztemp, args=(Om_m, Om_L))[0] / H0 #Mpc
			dAtemp = dLtemp / ((1 + ztemp) * (1 + ztemp)) #Mpc
	
			Ltemp = calcL(LOIIItemp)	
	
			#outf.write('z ' + str(ztemp) + '\n')
			#outf.write('L ' + str(Ltemp) + '\n')
			#outf.write('dL ' + str(dLtemp) + '\n')
			#outf.write('dist0 ' + str(dAtemp * 1000 * disttemp * math.pi / 180.) + '\n')

			dL.append(dLtemp)
			dA.append(dAtemp)
			z.append(ztemp)
			LOIII.append(LOIIItemp)
			L.append(Ltemp)
			dists.append(dAtemp * 1000 * disttemp * math.pi / 180.)

			mult = dAtemp * 1000 * lentemp * math.pi / 180.
			mult = mult / (ppix[i] - 1)

			for line in rawf:
				line = line.strip().split()
				line[0] = (float(line[0]) - 1) * mult
				line[1] = float(line[1]) / mval[i]
				outf.write(str(line[0]) + ' ' + str(line[1]) + '\n')

			rawf.close()
			outf.close()
	
		elif f.find('rad') >= 0:
			ztemp = float(rawf.readline().strip().split()[1])
			LOIIItemp = float(rawf.readline().strip().split()[1])
			lentemp = float(rawf.readline().strip().split()[1])

			Ltemp = calcL(LOIIItemp)

			#outf.write('z ' + str(ztemp) + '\n')
			#outf.write('L ' + str(Ltemp) + '\n')
			#outf.write('dL ' + str(dL[i]) + '\n')


			mult = dAtemp * 1000 * lentemp * math.pi / 180. 
			mult = mult/(rpix[i] - 1)

			for line in rawf:
				line = line.strip().split()
				line[0] = (float(line[0]) - 1) * mult
				line[1] = float(line[1]) / mval[i]
				outf.write(str(line[0]) + ' ' + str(line[1]) + '\n')

			i = i + 1
			rawf.close()
			outf.close()

		else:
			print 'This is crazy! Debug!'

	summary.write('\n z   \t')
	for j in range(len(z)):
		summary.write(str(z[j]).ljust(11) + '\t')
	summary.write('\n dL  \t')
	for j in range(len(dL)):
		summary.write(str(dL[j]).ljust(11) + '\t')
	summary.write('\n dA  \t')
	for j in range(len(dA)):
		summary.write(str(dA[j]).ljust(11) + '\t')
	summary.write('\n LOIII\t')
	for j in range(len(LOIII)):
		summary.write(str(LOIII[j]).ljust(11) + '\t')
	summary.write('\n L   \t')
	for j in range(len(L)):
		summary.write(str(L[j]) + '\t')
	summary.write('\n dists\t')
	for j in range(len(dists)):
		summary.write(str(dists[j]) + '\t')



main()

