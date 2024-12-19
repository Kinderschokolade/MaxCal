import numpy as np
import math
import pprint
from scipy.sparse.linalg.eigen.arpack import eigs

def J(x,y):
	if (x-y)==1: 
		return 1.
	elif (x-y)==-1: 
		return -1.
	else: 
		return 0

U = [4,1,1]
N=3

beta = 1./6.665 #one T profile, it should fix beta and nu, then vary U1,U2 to check
nu = 0#0.33 #constant for general motor

mu = -0.065#0.24 #const for single motor

pall = np.zeros((10,N))	
W = np.zeros((N,N))
P = np.zeros(N)

for k in range(0,10):
	U[0] = 0.5+ 0.5 *k	
	
	for i in range(N):
		for j in range(i-1,i+2):
			if j<0 : 
				j = 0
			if j == N : 
				j = N-1
	
			W[i][j] = math.exp(-beta * (U[i]+U[j])/2. - nu /N * J(i,j) )	
	
	W[0][N-1] =  math.exp(-beta * (U[0]+U[N-1])/2. - nu /N )	
	W[N-1][0] =  math.exp(-beta * (U[N-1]+U[0])/2. + nu /N )	
	
	W[N-1][0] = W[N-1][0] * math.exp(mu)
	W[N-1][1] = W[N-1][1] * math.exp(mu)
#	W[0][N-1] = W[0][N-1] * math.exp(-mu)
#	W[1][N-1] = W[1][N-1] * math.exp(-mu)

	Evl, Evecl = eigs(W, k=1, which='LR')
	Evr, Evecr = eigs(W.transpose(), k=1, which='LR')

	p = Evecl * Evecr
	p = p / p.sum()

	pall[k][0] = p[0][0].real
	pall[k][1] = p[1][0].real
	pall[k][2] = p[2][0].real


	
#pall = np.asarray(pall)
print(pall)

np.savetxt("/data/isilon/bause/single_particle/MaxCal_P_01.dat",pall)
	
