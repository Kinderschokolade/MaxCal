import numpy as np
import os
import math
import argparse

def H(x, kc, shift):
	return(1./(1.+np.exp(-2*kc*(x-shift))));



def dH(x, kc,shift):
	return(-2.*np.exp(-2.*kc*(-shift+x))*kc)/((1.+np.exp(-2*kc*(-shift+x)))*(1.+np.exp(-2*kc*(-shift+x))));


def F_tanh(x, kc, U, xc):
	lvl = np.size(U)
	f = (dH(x,kc,xc[0]) - dH(x,kc,xc[1])+ dH(x,kc,xc[lvl]))* U[0] +(dH(x,kc,xc[lvl-1]) - dH(x,kc,xc[lvl]) -dH(x,kc,xc[0]))* U[lvl-1]

	for i in range(1,lvl-1):
		f += (dH(x,kc,xc[i]) - dH(x,kc,xc[i+1]))* U[i]
	return f

def get_U(x, kc, U, xc):
	E = 0.
	lvl = len(U);

	E = (H(x,kc,xc[0]) - H(x,kc,xc[1])+ H(x,kc,xc[lvl]))* U[0] +(H(x,kc,xc[lvl-1]) - H(x,kc,xc[lvl]) - H(x,kc,xc[0]))* U[lvl-1] ;
	
	for j in range(1,lvl-1):
		E += (H(x,kc,xc[j]) - H(x,kc,xc[j+1]))* U[j] ;
		
	return E

def get_pos(x,ms):
	pos = math.floor(x*ms)
	if(pos == ms):
		 pos -= 1 # rare case
	return pos

parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
parser.add_argument('-n', help='start', default=7)
parser.add_argument('-m', help='end', default=17)
parser.add_argument('-ms', help='end', default=30)
args = parser.parse_args();
ident = args.o
n = int(args.n)
m = int(args.m)
ms = int(args.ms)




print(n,m)
path = "/data/pckr194/bause/single_particle/trajectory_" +str(ident) +".dat"
series = open(path, "r")
tra_len = int(os.popen("wc -l "+path).readline().split()[0])

U = np.zeros(6)
xc = np.zeros(7)

for line in open("/data/isilon/bause/single_particle/param_"+str(ident)+".dat","r"):
	cont = line.split()
	if(cont[0] == 'dT'):
		dT = float(cont[1])
	if(cont[0] == 'gamma'):
		gamma = float(cont[1])
	if(cont[0] == 'T0'):
		Temperature = float(cont[1])
	if(cont[0] == 'extf'):
		force = float(cont[1])
	if(cont[0] == 'U0'):
		U[0] = float(cont[1])
	if(cont[0] == 'U1'):
		U[1] = float(cont[1])
	if(cont[0] == 'U2'):
		U[2] = float(cont[1])
	if(cont[0] == 'U3'):
		U[3] = float(cont[1])
	if(cont[0] == 'U4'):
		U[4] = float(cont[1])
	if(cont[0] == 'U5'):
		U[5] = float(cont[1])
	if(cont[0] == 'xc0'):
		xc[0] = float(cont[1])
	if(cont[0] == 'xc1'):
		xc[1] = float(cont[1])
	if(cont[0] == 'xc2'):
		xc[2] = float(cont[1])
	if(cont[0] == 'xc3'):
		xc[3] = float(cont[1])
	if(cont[0] == 'xc4'):
		xc[4] = float(cont[1])
	if(cont[0] == 'xc5'):
		xc[5] = float(cont[1])
	if(cont[0] == 'xc6'):
		xc[6] = float(cont[1])
	if(cont[0] == 'k'):
		kc = float(cont[1])

stuff =series.readline()
split = stuff.split()
x = float(split[0])
pos = get_pos(x,ms)
pos_old = pos
Sprod = 0.
Sprod2 = 0.
SprodA = 0.

sqdT = np.sqrt(dT)
sigma = np.sqrt(2*Temperature/gamma)



traj = [] 
pos_l = []
Sprod_l = []
Sprod2_l = []
SprodA_l = []

lag =50
count =0
while (count < lag):
	stuff = series.readline()
	if stuff.strip():
		split = stuff.split()
		xold = x
		x = float(split[0])
		ran = float(split[2])
		F= F_tanh(x,kc,U,xc) + force
		pos_old = pos
		pos = get_pos(x,ms)
		traj.append(x)
		pos_l.append(pos)
		Sprod_l.append(2* F /(gamma*sigma) * sqdT *ran)
		Sprod2_l.append(2*F*F/(gamma*gamma*sigma*sigma) *dT)
		SprodA_l.append(F/Temperature * (xold-x) )
		Sprod += Sprod_l[count]
		Sprod2 += Sprod2_l[count]
		SprodA += SprodA_l[count]
		count = count +1

	
start = m
end = n

count = 0
ccount = 0
path1 = "/data/pckr194/bause/traj/tr1_"
path2 = "/data/pckr194/bause/traj/tr2_"

#tra_len = 100000
potential = np.loadtxt("/data/isilon/bause/single_particle/potential_ms_"+str(ident)+".dat")

print ( m /ms , n/ms)
Sprod_th = ( -force / ms * (m-n) + get_U(m/ms,kc,U,xc) - get_U(n/ms,kc,U,xc)) / Temperature
print("theory" ,Sprod_th)

Sprodf_data= []
Sprodb_data= []
SprodAb_data= []
SprodAf_data= []


#tlen = 0

ran_l =[0.]*lag


for i in range(tra_len):
	stuff = series.readline()
	if stuff.strip():	
		split = stuff.split()
		xold = x
		x = float(split[0])
		ran = float(split[2])
		F= F_tanh(x,kc,U,xc) + force
		pos_old = pos
		pos = get_pos(x,ms)
		traj.append(x)
		ran_l.append(ran)
		if (pos == m and start == n):
			start = m
			#SprodAf_data.append(Sprod+Sprod2)
#			Sprodg = Sprod+Sprod2
			length = len(traj)
			Sprod =0
			Sprod2 =0
#			Sprod_alt = 0
			count = count +1 
			for j in range(length):
				F= +F_tanh(traj[length-j-1],kc,U,xc)+ F_tanh(traj[j],kc,U,xc) - 2.*force 
				Sprod+= -F /(gamma*sigma) * sqdT *ran_l[j]
				Sprod2+= F*F/(2*gamma*gamma*sigma*sigma) *dT
			
			Sprodg = Sprod+ Sprod2				
			if (length < 81):	
				Sprodf_data.append(Sprodg)
#				print(count,Sprod_alt,Sprodg,Sprod,Sprod2,Sprod_th,len(traj))
		#		np.savetxt(path1+str(count)+".dat",np.asarray(traj))		
			traj = [x]
			ran_l = [ran]
			Sprodg = Sprod+ Sprod2				
			SprodAf_data.append(Sprodg)
#	
#			if((Sprod+Sprod2) > 1.*Sprod_th):
#				print(count,Sprod+Sprod2,Sprod,Sprod2,Sprod_th)
		if (pos == n and start == m):
			start = n
			#SprodAf_data.append(Sprod+Sprod2)
#			Sprodg = Sprod+Sprod2
			ccount = ccount +1;
			length = len(traj)
			for j in range(length):
				F= +F_tanh(traj[length-j-1],kc,U,xc) + F_tanh(traj[j],kc,U,xc) - 2.*force  # - force from inveresed trajectory
				Sprod+=  -F /(gamma*sigma) * sqdT *ran_l[j]
				Sprod2+= F*F/(gamma*gamma*sigma*sigma*2.) *dT

			Sprodg = Sprod+ Sprod2				
			traj = [x]
			ran_l = [ran]
			Sprodg = Sprod+Sprod2
			SprodAb_data.append(Sprodg)
			if (length < 81):
				Sprodb_data.append(Sprodg)
			#print(ccount,Sprodg,Sprod,Sprod2,Sprod_th)
			#np.savetxt(path2+str(ccount)+".dat",np.asarray(traj))		

#print("min: ", pos_min , Sprod_min, Sprod_th)
#print("max: ", pos_max , Sprod_max, Sprod_th)
#
#
bins = 100
ra = 4* abs(Sprod_th)
histf,edgesf = np.histogram(np.asarray(Sprodf_data),bins = bins,range=(-ra,ra))
histb,edgesb = np.histogram(np.asarray(Sprodb_data),bins = bins,range=(-ra,ra))
histAf,edgesAf = np.histogram(np.asarray(SprodAf_data),bins = bins,range=(-ra,ra))
histAb,edgesAb = np.histogram(np.asarray(SprodAb_data),bins = bins,range=(-ra,ra))
histff = np.zeros((bins,2))
histbb = np.zeros((bins,2))
histAbb = np.zeros((bins,2))
histAff = np.zeros((bins,2))
#
for i in range(bins):
	histff[i,0] = (edgesf[i] + edgesf[i+1]) / 2.
	histbb[i,0] = (edgesb[i] + edgesb[i+1]) / 2.
	histAff[i,0] = (edgesAf[i] + edgesAf[i+1]) / 2.
	histAbb[i,0] = (edgesAb[i] + edgesAb[i+1]) / 2.
histff[:,1] = histf
histbb[:,1] = histb
histAff[:,1] = histAf
histAbb[:,1] = histAb
#
##meanf = np.inner(histff[:,0],histff[:,1])/ sum(histf)
##meanb = np.dot(histbb[:,0],histbb[:,1])/ sum(histb)
#
##print(meanf, meanb)
#
np.savetxt("/data/isilon/bause/single_particle/Sprod/hist_"+str(m)+"_"+str(n)+"_"+str(ident)+".dat", histff)
np.savetxt("/data/isilon/bause/single_particle/Sprod/hist_"+str(n)+"_"+str(m)+"_"+str(ident)+".dat", histbb)
np.savetxt("/data/isilon/bause/single_particle/Sprod/histA_"+str(n)+"_"+str(m)+"_"+str(ident)+".dat", histAbb)
np.savetxt("/data/isilon/bause/single_particle/Sprod/histA_"+str(m)+"_"+str(n)+"_"+str(ident)+".dat", histAff)
#
#
