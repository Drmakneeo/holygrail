import math

from sympy import *

from astropy import units as u

from poliastro.bodies import Moon
from poliastro.twobody import Orbit
from poliastro import constants

from scipy.special import lpmn

C = [[0], [0]]
S = [[0], [0]]

with open('glgm2sh.tab') as fp:
	line = fp.readline()
	while line:
		val = line.split(',')
		m = int(val[0])
		n = int(val[1])
		if (n == 0):
			C.append([float(val[2])])
			S.append([float(val[3])])
		else:
			C[m].append(float(val[2]))
			S[m].append(float(val[3]))
		line = fp.readline()

P = [[0], [0]]
z = Symbol('z')
for m in range(2, 11):
	f = ((z**2) - 1)**m
	print(f)
	for i in range(0, m):
		f = f.diff(z)
		print(f)
	for n in range(0, m):
		l = ((-1)**n)/((math.factorial(m)*2**m))*((1-z**2)**(m/2))
		if (n == 0):
			s = f * l
			P.append([s])
		else:
			f = f.diff(z) * l
			print(s)
			P[m].append(s)

def gravity(r, theta, phi):
	V = constants.GM_moon / (constants.R_mean_moon**2)
#	z = math.cos(theta)
#	P = lpmn(70, 70, z)[0]
	x = 1
	for m in range(2, 11):
		for n in range(0, m):
			Pmn = lambdify(z, P[m][n])
			print(Pmn(math.cos(theta)))
			s1 = Pmn(math.cos(theta))
			s2 = s1 * (S[m][n])*math.sin(m*phi)
			s3 = s1 * (C[m][n])*math.cos(m*phi)
			s = (s2 + s3) * (constants.R_mean_moon / (r*u.m))**n
			x += s
	return V * x

def main():
	print(gravity(2e7, 2, 4))

if __name__ == '__main__':
    main()