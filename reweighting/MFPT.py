import numpy as np
import argparse
import getpass
import os.path
from define_flow import define_flow_direction
from define_flow import read_cut
from define_flow import read_minpos
import math
from msmtools.analysis import mfpt
from fpt_ana import  first_passage



def mfpt_markov(L, minpos, rancut, ms,lvl):
        mfpt_markov = np.zeros((lvl,lvl))
        for k in range(lvl):
                target = list(range(minpos[k] - rancut , minpos[k] + rancut))
                for t in range(lvl):
                        origin = list(range(minpos[t] - rancut , minpos[t] + rancut))
                        #print(origin, target)
                        if(origin[-1] >= ms or target[-1] >= ms):
                                print ("origin/target for mfpt out of range!")
                        else:
                                mfpt_markov[t,k] = mfpt(L,target = target, origin = origin)
        return mfpt_markov


def reweight(T,r, gamma, ms):
        W = [ [ T[i,j] * math.exp(gamma*r[i,j]) for j in range(ms)] for i in range(ms)]
        W = np.asarray(W)

        Ev,Evecl = np.linalg.eig(np.transpose(W))
        idx = Ev.argsort()[::-1]
        Ev = Ev[idx]  #order eigenvectors by size
        Evecl = np.transpose(Evecl)
        Evecl = Evecl[idx,:]
        Ev,Evecr = np.linalg.eig(W)
        idx = Ev.argsort()[::-1]
        Ev = Ev[idx]  #order eigenvectors by size
        Evecr = np.transpose(Evecr)
        Evecr = Evecr[idx,:]
        p  = Evecl[0,:] * Evecr[0,:]
        p = np.absolute(np.asarray(p)/sum(p))
        k = [[ Evecr[0,j] / (Evecr[0,i] * Ev[0]) *W[i,j] for j in range(ms)] for i in range(ms)]
        k= np.asarray(k)

        for i in range(ms):
                k[i] = k[i] / sum(k[i,:])

        return np.real(k),np.real(p)


def integrate_MFPT(FPT,length):
	res = 0.
	for i in range(length):
		res+= i*FPT[i]
	return res

def calc_av(F, T, p):
	ms = np.shape(T)[0]
	return np.real( sum( sum( p[i] * T[i,j] *F[i,j] for j in range(ms)) for i in range(ms)  ))

def analyse_MSM(k):
	Ev,Evec = np.linalg.eig(np.transpose(k))
	idx = Ev.argsort()[::-1]
	Ev = Ev[idx]  #order eigenvectors by size
	Evec = np.transpose(Evec)
	Evec = Evec[idx,:]
	p = Evec[0] / sum(Evec[0])
	return Ev,p,Evec



parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident start')
parser.add_argument('-b', help='number of sets, count up x (kxkk)')
parser.add_argument('-c', help='number of sets, count up x (kxkk)')
parser.add_argument('-l', help='lagtime chosen')
args = parser.parse_args();
ident = int(args.o)
cmax = int(args.c)  
bmax = int(args.b)
l = int(args.l)

user = getpass.getuser()

if (user=="mbause"):
        path0 = "/u/mbause/data/"
        path1 = "/u/mbause/data/"
else:
        path0 = "/data/isilon/bause/"
        path1 = "/data/pckr194/bause/"

ms,cut = read_cut(path0,ident)
ms, minpos = read_minpos(path0,ident)
F = define_flow_direction(ms,cut)

MFPT = []
J = []
MFPT_er = []
MFPT_sig = []
J_er = []
rancut =1
lvl =3 

for i in range(bmax):
	av = np.zeros(lvl*lvl)
	av2= np.zeros(lvl*lvl)
	avJ = 0
	av2J =  0
	count = 0
	for j in range(cmax):
		ID = str(ident+ i * 100 + j*10)
		filepath= path0+"single_particle/"+str(l)+"/MFPT/markovana_"+ID+".dat"
		if (os.path.isfile(filepath)):
			data = np.loadtxt(filepath)
			T = np.loadtxt(path0+"single_particle/"+str(l)+"/T/T_"+ID+".dat")
			EV,p,Evec = analyse_MSM(T)
			#res = integrate_MFPT(data[:,1],300)
			res = mfpt_markov(T, minpos, rancut, ms,lvl).reshape(lvl*lvl)
			av += res
			av2 += res*res
			count +=1
			flux = calc_av(F,T,p)
			avJ += flux
			av2J += flux*flux
	av = av / count
	av2 = av2/ count
	MFPT.append(av)
	mfpt_temp = av2-av*av
	MFPT_sig.append(mfpt_temp)
	MFPT_er.append(mfpt_temp / np.sqrt(count))
	avJ = avJ / count
	av2J = av2J /count
	J.append(avJ)
	J_er.append((av2J-avJ*avJ) / np.sqrt(count))


J = np.transpose([J])
J_er = np.transpose([J_er])
#MFPT = np.transpose(MFPT)
#MFPT_er = np.transpose(MFPT_er)
#MFPT_sig = np.transpose(MFPT_sig)
MFPT = np.asarray(MFPT)
MFPT_er = np.asarray(MFPT_er)
MFPT_sig = np.asarray(MFPT_sig)
print(J.shape, J_er.shape, MFPT.shape, MFPT_er.shape)
out = np.concatenate((J,MFPT,MFPT_sig,MFPT_er),axis =1 )
		
np.savetxt(path0+"single_particle/"+str(l)+"/MFPT_"+str(ident)+".dat",out)



