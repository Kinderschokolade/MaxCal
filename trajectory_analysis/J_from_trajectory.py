import numpy as np
import  argparse
import os
from define_flow import read_minpos
from define_flow import read_cut


parser = argparse.ArgumentParser()
parser.add_argument('-o', help='ident')
args = parser.parse_args();
ident = args.o


path = "/data/isilon/bause/"
ms, minpos = read_minpos(path,ident)
ms, cut = read_cut(path,ident)

L = np.zeros((ms,ms))
lvl = len(minpos)
mfpt_ar = np.zeros((lvl,lvl))


path = "/data/pckr194/bause/single_particle/trajectory_" +str(ident) +".dat"
series = open(path, "r")
tra_len = int(os.popen("wc -l "+path).readline().split()[0])


stuff =series.readline()
split = stuff.split()
x = float(split[0])

num_block = 10
winding = 0;
dT =  0.00001
s =0
trav = np.zeros(num_block)

tra_len_bl = int(tra_len / num_block)

for j in range(num_block):
	s = 0;
	for i in range(tra_len_bl):
		stuff = series.readline()
		if stuff.strip():
			x_old = x
			split = stuff.split()
			x = float(split[0])
			diff = abs(x - x_old)
			s +=1	

			if (diff > 0.5):
				if (x_old > x):
					winding +=1
					trav[j]+= x_old + x -1.
				else:
					winding -=1
					trav[j]+= -x_old-x + 1.
			else:
				trav[j] += x - x_old

	trav[j]= trav[j] / (dT*s)	
			
J = sum(trav) / num_block
trav2 =  sum(trav[i] * trav[i] /num_block for i in range(num_block) )
Jsigma = trav2 - J*J 

print(ident, J, Jsigma)





